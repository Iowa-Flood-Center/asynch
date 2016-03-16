/* jshint node: true, esnext: true */

"use strict";

var fs = require('fs-extra');
var path = require('path');
var cp = require('child_process');
var os = require('os');

var username = require('username');
var parseArgs = require('minimist');
var debug = require('debug')('forecaster');
var pgp = require('pg-promise')({
  // Initialization Options
});

var templates = require('./templates.js');

// Parse the command line arguments
var argv = parseArgs(process.argv.slice(2), {
  string: ['model', 'output'],
  alias: {
    o: 'output'
  },
  default: {
    output: 'out'
  }
});

// Read the config file
try {
  var cfg = JSON.parse(fs.readFileSync('.forecaster', 'utf8'));
} catch (ex) {
  debug('Missing or invalid config file found, generating one');
  debug(ex.message);
  
  // Get current time rounded to the nearest hour  
  var now = Math.floor(Date.now() / 1000);
  now = now - (now % 3600);
  
  var cfg = {
    username: username.sync(),
    obsTime: now,
    qpfTime: now
  };
  fs.writeFileSync('.forecaster', JSON.stringify(cfg), 'utf8');
}

if (cfg && cfg.jobId) {
  try {
    var stat = cp.execSync('qstat -j ' + cfg.jobId).toString();
    debug(stat);
  } catch(ex) {
    debug(ex);
  }
}

// Set the ouput directory
var outputDir = argv.output;

// Create output dir if not exists
fs.mkdirsSync(outputDir);

function renderAll(context) {
  // Render all templates
  Object.keys(templates).forEach(function (key, index) {

    var file = path.join(outputDir, prefix + key);
    debug('Rendering ' + file);

    var content = templates[key](context);
    fs.writeFileSync(file, content);
  });
}

function render(template, out, context) {
  var file = path.join(outputDir, out);
  debug('Rendering ' + file);

  var content = template(context);
  fs.writeFileSync(file, content);
}

function generateStormfile(filePath, rows) {  
  //Maps the link
  var links = {};

  var length = rows.length;
  for (var i = 0 ; i < length ; i++) {
    if (typeof links[rows[i].link_id] === 'undefined') {
      links[rows[i].link_id] = [{
        time: rows.length,
        value: rows.rain_intens
      }];
    } else {
      links[rows[i].link_id].push({
        time: rows.length,
        value: rows.rain_intens
      });
    }    
  }
  
  // Generate the storm files
  var outFile = fs.createWriteStream(filePath);
  outFile.write(Object.keys(links).length +  os.EOL);
  
  Object.keys(links).forEach(function(key, index) {
    var link = links[key];
    outFile.write(key + ' ' + link.length + os.EOL);
    link.forEach(function(change) {
      outFile.write(change.time + ' ' + change.value + os.EOL);
    });
  });
                             
  outFile.end();           
}

var obsConnString = 'postgres://***REMOVED***:***REMOVED***@s-iihr51.iihr.uiowa.edu/model_ifc';
var forecastConnString = 'postgres://***REMOVED***:***REMOVED***@s-iihr51.iihr.uiowa.edu/h3r_qpf';

// Creating a new database instance from the connection details
var dbObs = pgp(obsConnString),
  dbQpf = pgp(forecastConnString);

// Query the DB to get the latest OBS and QPF timestamp
Promise.all([
  dbObs.one('SELECT unix_time FROM rain_maps5_index ORDER BY unix_time DESC LIMIT 1'),
  dbQpf.one('SELECT time_utc FROM hrrr_index_ldm ORDER BY time_utc DESC LIMIT 1')
])
  .then(function (rows) {
    var obsTime = rows[0].unix_time;
    var qpfTime = rows[1].time_utc;
  
    debug('lastest rainfall OBS timestamp ' + obsTime);
    debug('lastest rainfall QPF timestamp ' + qpfTime);

    // If we have new obs or qpf
    if ((cfg.obsTime < obsTime) || (cfg.qpfTime < qpfTime)) {
      var obsContext = {
          begin: cfg.obsTime,
          end: obsTime,
          duration: obsTime - cfg.obsTime,
          user: cfg.username,
          rainFile: 'forcing_rain_obs.str',
          stateFile: 'obs.rec'
        },
        qpfContext = {
          begin: obsTime,
          end: qpfTime,
          duration: qpfTime - obsTime,
          user: cfg.username,
          rainFile: 'forcing_rain_qpf.str',
          stateFile: 'qpf.rec'
        };
      
      // Render a new set of config files
      render(templates['gbl'], 'obs.gbl', obsContext);
      render(templates['gbl'], 'qpf.gbl', qpfContext);
      render(templates['job'],  'run.job', {
          globalFiles: ['obs.gbl', 'qpf.gbl'],
          workingDir: path.resolve(outputDir)	
        });

      debug('select OBS from ' + obsContext.begin + ' to ' + obsContext.end);
      debug('select QPF from ' + qpfContext.begin + ' to ' + qpfContext.end);
      
      // Generate forcing files
      Promise.all([
        dbObs.any('SELECT unix_time, rain_intens, link_id FROM link_rain5 ' +
                  'WHERE unix_time >= $1 AND unix_time < $2 ' +
                  'ORDER BY link_id, unix_time',
                  [obsContext.begin, obsContext.end])
          .then(function(rows) {
            debug('got ' + rows.length + ' OBS rainfall rows');
            generateStormfile(path.join(outputDir, obsContext.rainFile), rows);
            //dbObs.done();
          }),
        dbQpf.any('SELECT unix_time, rain_intens, link_id FROM link_rain ' +
                   'WHERE unix_time >= $1 AND unix_time < $2 ' +                        
                   'ORDER BY link_id, unix_time',
                   [qpfContext.begin, qpfContext.end])
          .then(function(rows) {
            debug('got ' + rows.length + ' QPF rainfall rows');
            generateStormfile(path.join(outputDir, qpfContext.rainFile), rows);
            //dbQpf.done();
          })
        ]).then(function() {
        
          // Run the simulations
          debug('Run the simulations');
          try {
            var re = /Your job (\d*)/;

            function getJobId(sub) {
              var m = re.exec(sub);
              if (m !== null) return parseInt(m[1]);
            }

            // Queue the job
            var sub = cp.execSync('qsub ' + path.join(outputDir, 'run.job')).toString();
            debug(sub);
 
            // Update config file for next run
            cfg.jobId = getJobId(sub);
            cfg.obsTime = obsTime;
            cfg.qpfTime = qpfTime;

            fs.writeFileSync('.forecaster', JSON.stringify(cfg), 'utf8');
          } catch (ex) {
            console.error(ex.message);
          }
        }).catch(function (err) {
          return console.error('error running query', err);
        });
    }

  })
  .catch(function (err) {
    return console.error('error running query', err);
  });

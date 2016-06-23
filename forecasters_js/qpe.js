/* jshint node: true, esnext: true */

"use strict";

var fs = require('fs-extra');
var path = require('path');
var os = require('os');

var username = require('username');
var parseArgs = require('minimist');
var debug = require('debug')('forecaster');
var pgp = require('pg-promise-strict');

var config = require('./config.js');
var templates = require('./templates.js');
var sge = require('./sge.js');
var stormfile = require('./stormfile.js');
//var scenarii = require('./scenarii.js');

// Parse the command line arguments
var argv = parseArgs(process.argv.slice(2), {
  string: ['output'],
  boolean: 'qsub',
  alias: {
    o: 'output'
  },
  default: {
    output: 'out'
  }
});

// Set the ouput directory
var outputDir = argv.output;

// Get latest file matching a given pattern
function getLatest(re) {
  return fs.readdirSync(outputDir)
    .filter(function (filename) {
      var stat = fs.statSync(path.join(outputDir, filename));
      return (stat.isFile() && re.test(filename));
    })
    .map(function (basename) {
      var m = re.exec(basename);
      if (m !== null) {
        return {
          filename: m[0],
          time: parseInt(m[1])
        };
      }
      return null;
    })
    .reduce(function(max, x) {
      return x.time > max.time ? x : max;
    });
}

// Get the latest available state
function getLatestState() {
  return getLatest(/^state_(\d+).(rec|h5)$/);
}

// Get the latest forcing
function getLatestForcing() {
  return getLatest(/^forcing_rain_qpe_(\d+).str$/);
}

// Create output dir if not exists
fs.mkdirsSync(outputDir);

//Render a template
function render(template, out, context) {
  var file = path.join(outputDir, out);
  debug('Rendering ' + file);

  var content = template(context);
  fs.writeFileSync(file, content);
}

//Check whether a simulation is already running
sge.qstat('IFC_QPE')
.then(function (status) {
  debug(status);
  if (status && status.state.indexOf('r') !== -1) {
    throw 'Simulation already running';
  } else if (status && status.state.indexOf('w') !== -1) {
    debug('Simulation is waiting, put on hold while job is updated');
    return sge.qhold('IFC_QPE');
  } else {
    debug('No job found, submit a new one');
    render(templates.job, 'qpe.job', {
      name: 'IFC_QPE',
      globalFile: 'qpe.gbl',
      workingDir: path.resolve(outputDir)
    });

    return sge.qsub(path.join(outputDir, 'qpe.job'));
  }
})
.then(function () {
  // Creating a new database instance from the connection details
  return pgp.connect(config.qpeConn).then(function (client) {
    // Query the DB to get the latest obs timestamp
    return client.query('SELECT unix_time FROM rain_maps5_index ORDER BY unix_time DESC LIMIT 1');
  });
})
.then(function (query){
  return query.fetchUniqueValue();
}).then(function (result) {
  var latestState = getLatestState(),
      latestForcing = getLatestForcing();
  
  var currentTime = latestState.time,
    forcingTime = latestForcing.time || currentTime,
    latestTime = result.value,
    user = username.sync();
  
  var iniStateMode;
  switch(path.extname(latestState.filename)) {
    case '.rec':
      iniStateMode = 2;
      break;
    case '.h5':
      iniStateMode = 4;
      break;
    default:
      console.error('Unkown initial state file extension');
      process.exit(1);
  }

  debug('current timestamp ' + currentTime);
  debug('last rainfall QPE timestamp ' + forcingTime);
  debug('lastest rainfall QPE timestamp ' + latestTime);

  // If we have new obs ready
  if (forcingTime >= latestTime) {
    debug ('Simulation is already up to date');
    return sge.qrls('IFC_QPE');
  }

  var context = {
    user: user,
    begin: currentTime,
    end: latestTime ,
    duration: (latestTime - currentTime) / 60,
    iniStateMode: iniStateMode,
    iniStateFile: latestState.filename,
    endStateFile: 'state_' + latestTime + '.h5',
    rainFileType: 1,
    rainFile: latestForcing.filename,
    outHydrographsFile: {
      file: 'hydrographs.dat'
    },
    outPeaksFile: {
      file: 'peaks.pea'
    }
  };

  // Render a new set of config files
  render(templates.gbl, 'qpe.gbl', context);

  // Get the QPE data and generate the stormfile
  var links = {};
  return result.client.query(config.qpeQuery, [context.begin])
    .onRow(stormfile.mapLink(context.begin, 300, links))
  .then(function (result) {
    debug('got ' + result.rowCount + ' QPE rainfall rows');
    result.client.done();

    stormfile.generate(path.join(outputDir, context.rainFile), links);

    // Run the simulations
    debug('Release the simulation');
    return sge.qrls('IFC_QPE');
  });
})
.catch(function (err) {
  return debug(err);
}).then(function(){
  process.exit();
});
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
    }, {time: null});
}

// Get the latest available state
function getLatestState() {
  return getLatest(/^state_ifc_(\d+).(rec|h5)$/);
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
sge.qstat('WHATIF_NORAIN')
.then(function (status) {
  debug(status);
  if (status && status.state.indexOf('r') !== -1) {
    throw 'Simulation already running';
  } else if (status && status.state.indexOf('w') !== -1) {
    debug('Simulation is waiting, put on hold while job is updated');
    return sge.qhold('WHATIF_NORAIN');
  } else {
    debug('No job found, submit a new one');
    render(templates.job, 'whatif_norain.job', {
      name: 'WHATIF_NORAIN',
      globalFile: 'whatif_norain.gbl',
      workingDir: path.resolve(outputDir)
    });

    return sge.qsub(path.join(outputDir, 'whatif_norain.job'), ['QPE_IFC']);
  }
})
.then(function () {
  var latestState = getLatestState();
  
  var currentTime = latestState.time,
    endTime = currentTime + 14400 * 60,
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

  var context = {
    user: user,
    begin: currentTime,
    end: endTime,
    duration: 14400,
    rainFileType: 0,
    iniStateMode: iniStateMode,
    iniStateFile: latestState.filename,
    endStateFile: 'forecast_norain_' + endTime + '.h5',
    outHydrographsDb: {
      file: 'save_hydrograph.dbc',
      table: 'whatif_norain.hydrographs'
    },
    outPeaksDb: {
      file: 'save_peaks.dbc',
      table: 'whatif_norain.peaks'
    }
  };

  // Render a new set of config files
  render(templates.gbl, 'whatif_norain.gbl', context);

  // Run the simulations
  debug('Release the simulation');
  sge.qrls('IFC_NORAIN');

  //Update the run start time
  return pgp.connect(config.outConn).then(function (client) {
    debug('updating run start time ' + context.begin);
    return client.query('UPDATE runs SET start_time = to_timestamp($1) WHERE run LIKE $2', [context.begin, 'whatif_norain']);
  })
  .then(function (query){
    return query.fetchOneRowIfExists();
  })
  .then(function (result) {
    result.client.done();
    debug('done.');
  })
  .catch(function (err) {
    debug(err);
  });
})
.catch(function (err) {
  return debug(err);
}).then(function(){
  process.exit();
});
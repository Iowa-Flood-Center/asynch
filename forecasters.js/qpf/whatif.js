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
  var dates = fs.readdirSync(outputDir)
    .filter(function (filename) {
      var stat = fs.statSync(path.join(outputDir, filename));
      return (stat.isFile() && re.test(filename));
    })
    .map(function (basename) {
      var m = re.exec(basename);
      if (m !== null) {
        return parseInt(m[1]);
      }
      return null;
    });

  return Math.max.apply(Math, dates);
}

// Get the latest available state
function getLatestState() {
  return getLatest(/^state_(\d+).rec$/);
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
sge.qstat('IFC_WHATIF_2IN24H')
.then(function (status) {
  debug(status);
  if (status && status.state.indexOf('r') !== -1) {
    throw 'Simulation already running';
  } else if (status && status.state.indexOf('w') !== -1) {
    debug('Simulation is waiting, put on hold while job is updated');
    return sge.qhold('IFC_WHATIF_2IN24H');
  } else {
    debug('No job found, submit a new one');
    render(templates.job, 'whatif_2in24h.job', {
      name: 'IFC_WHATIF_2IN24H',
      globalFile: 'whatif_2in24h.gbl',
      workingDir: path.resolve(outputDir)
    });

    return sge.qsub(path.join(outputDir, 'whatif_2in24h.job'), ['IFC_QPE']);
  }
})
.then(function () {
  var currentTime = getLatestState(),
    user = username.sync();

  debug('current timestamp ' + currentTime);

  var context = {
    user: user,
    begin: currentTime,
    end: currentTime + 14400 * 60,
    duration: 14400,
    rainFileType: 4,
    rainFile: 'forcing_rain_2inches24hour.ustr',
    iniStateFile: 'state_' + currentTime + '.rec',
    endStateFile: 'forecast_2in24h_.rec',
    outHydrographsDb: {
      file: 'save_hydrograph.dbc',
      table: 'whatif_2in24h.hydrographs'
    },
    outPeaksDb: {
      file: 'save_peaks.dbc',
      table: 'whatif_2in24h.peaks'
    }
  };

  // Render a new set of config files
  render(templates.gbl, 'whatif_2in24h.gbl', context);

  // Run the simulations
  debug('Release the simulation');
  return sge.qrls('IFC_WHATIF_2IN24H');
})
.catch(function (err) {
  return debug(err);
}).then(function(){
  process.exit();
});
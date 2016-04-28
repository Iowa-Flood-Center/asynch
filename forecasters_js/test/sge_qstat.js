var debug = require('debug')('sge');

var sge = require('../sge.js');

sge.qstat('IFC_QPE')
.then(function (status) {
  debug(status);
  if (status && status.state.indexOf('r') !== -1) {
    throw 'Simulation already running';
  } else if (status && status.state.indexOf('w') !== -1) {
    debug('Simulation is waiting');
  } else {
    return status;
  }
})
.catch(function (err) {
  return console.error(err);
});

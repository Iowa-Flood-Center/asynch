var debug = require('debug')('sge');

var sge = require('../sge.js');

sge.qrls('IFC_QPE')
.catch(function (err) {
  return debug(err);
});

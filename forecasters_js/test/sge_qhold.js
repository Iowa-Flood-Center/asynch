var debug = require('debug')('sge');

var sge = require('../sge.js');

sge.qhold('IFC_QPE')
.catch(function (err) {
  return debug(err);
});

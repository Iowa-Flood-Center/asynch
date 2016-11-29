/* jshint node: true, esnext: true */

"use strict";

var fs = require('fs');
var os = require('os');

module.exports = {

  // For each rows, add value to the link rain data
  mapLink: function (begin, dt, links) {

    return function (row) {
      if (typeof links[row.link_id] === 'undefined') {
        links[row.link_id] = {};
        links[row.link_id][0] = 0.0;
      }

      if (row.unix_time !== null) {
        // Timestamps are at the end of the time step
        var t = (Math.max(row.unix_time - dt, begin) - begin) / 60;
        links[row.link_id][t] = row.rain_intens;
        links[row.link_id][t + Math.ceil(t / dt) * dt / 60] = 0.0;
      }
    };
  },

  // Generate the storm files
  generate: function (filePath, links) {
    // Generate the storm files
    var outFile = fs.createWriteStream(filePath);
    outFile.write(Object.keys(links).length + os.EOL);

    Object.keys(links).forEach(function (key, index) {
      var link = links[key];
      outFile.write(key + ' ' + Object.keys(link).length + os.EOL);
      Object.keys(link).forEach(function (time, index) {
        outFile.write(time + ' ' + link[time] + os.EOL);
      });
    });

    outFile.end();
  }

};
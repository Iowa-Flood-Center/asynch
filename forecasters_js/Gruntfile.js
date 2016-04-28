/*jslint node: true */

module.exports = function (grunt) {
  'use strict';

  // Time how long tasks take. Can help when optimizing build times
  require('time-grunt')(grunt);
  require('load-grunt-tasks')(grunt);

    // Define the configuration for all the tasks
  grunt.initConfig({

        // Precompile Hogan templates to javascript
    hogan: {
      options: {
        binderName: 'nodejs',
        nameFunc: function (filename) {
          var path = require('path');
          var basename = path
            .basename(filename, '.mustache');

          // to mixedCase
          return basename.replace(/-([a-z])/g, function (g) { return g[1].toUpperCase(); });
        }
      },
      app: {
        src: 'templates/*.*',
        dest: 'templates.js'
      }
    },

  });

	grunt.registerTask('default', [
    'hogan'
  ]);
};

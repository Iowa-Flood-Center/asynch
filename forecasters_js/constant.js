module.exports = {
  obsConn: 'postgres://automated_solver:C5.pfest0@s-iihr51.iihr.uiowa.edu/model_ifc',
  qpfConn: 'postgres://automated_solver:C5.pfest0@s-iihr51.iihr.uiowa.edu/h3r_qpf',
  obsQuery: `
SELECT unix_time, rain_intens, link_id
FROM
(SELECT link_id FROM materialized_env_master_km) links LEFT JOIN
(SELECT * FROM link_rain5 WHERE unix_time > $1 AND rain_intens > 0.0) rain USING (link_id)
ORDER BY link_id, unix_time`,
  qpfQuery: `
SELECT unix_time, rain_intens, link_id
FROM
(SELECT link_id FROM materialized_env_master_km) links LEFT JOIN
(SELECT * FROM link_rain WHERE unix_time > $1 AND rain_intens > 0.0) rain USING (link_id)
ORDER BY link_id, unix_time`
};

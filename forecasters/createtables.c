//Creates tables and functions in the rm_model database.
//gcc createtables.c -o CREATETABLES -O3 -lpq
#include <stdio.h>
#include <stdlib.h>
#include <libpq-fe.h>
#include <string.h>

int CheckSQLError(PGresult* res);
void CheckConnConnection(PGconn* conn);

int main(int argc,char* argv[])
{
	int i,j,numtables,program,new_version;
	char query[16384];
	char M[32];
	PGresult *res;
	PGconn *conn;

	if(argc < 3)
	{
		printf("Need model name, program type (forecast, maps).\n");
		return 1;
	}

	if(strcmp(argv[2],"forecast") == 0)	program = 0;
	else if(strcmp(argv[2],"maps") == 0)	program = 1;
	else
	{
		printf("Bad program type %s.\n",argv[2]);
		return 1;
	}

	//Connect to database
	conn = PQconnectdb("dbname=blah host=blah port=blah user=blah password=blah"); new_version = 1;
	CheckConnConnection(conn);

	//Load up model names
	sprintf(M,argv[1]);
	numtables = 10;
	printf("Using model name %s.\nAssuming %i archive tables.\n\n",M,numtables);

	//Create archive tables
	printf("Creating tables...\n");

	if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS master_archive_hydroforecast_%s (link_id integer,time_utc timestamp with time zone,discharge double precision,forecast_time integer,baseflow double precision);",M);
	else		sprintf(query,"CREATE TABLE master_archive_hydroforecast_%s (link_id integer,time_utc timestamp with time zone,discharge double precision,forecast_time integer,baseflow double precision);",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	for(i=0;i<numtables;i++)
	{
		if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS archive_hydroforecast_%s_%i() INHERITS (master_archive_hydroforecast_%s);",M,i,M,M,i);
		else		sprintf(query,"CREATE TABLE archive_hydroforecast_%s_%i() INHERITS (master_archive_hydroforecast_%s);",M,i,M,M,i);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);

		sprintf(query,"CREATE INDEX idx_archive_hydroforecast_%s_%i_forecast_time_link_id ON archive_hydroforecast_%s_%i USING btree (forecast_time, link_id);",M,i,M,i);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);
	}

	if(program == 1)
	{
		if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS master_archive_peakflows_%s (link_id integer,peak_time integer,peak_discharge double precision,forecast_time integer,period integer);",M);
		else		sprintf(query,"CREATE TABLE master_archive_peakflows_%s (link_id integer,peak_time integer,peak_discharge double precision,forecast_time integer,period integer);",M);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);

		for(i=0;i<numtables;i++)
		{
			if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS archive_peakflows_%s_%i() INHERITS (master_archive_peakflows_%s);",M,i,M,M,i);
			else		sprintf(query,"CREATE TABLE archive_peakflows_%s_%i() INHERITS (master_archive_peakflows_%s);",M,i,M,M,i);
			res = PQexec(conn,query);
			CheckSQLError(res);
			PQclear(res);

			sprintf(query,"CREATE INDEX idx_archive_peakflows_%s_%i_forecast_time_link_id ON archive_peakflows_%s_%i USING btree (forecast_time, link_id);",M,i,M,i);
			res = PQexec(conn,query);
			CheckSQLError(res);
			PQclear(res);
		}

		if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS master_archive_maps_%s (forecast_time integer,link_id integer,q double precision,S double precision,s_p double precision,s_l double precision,s_s double precision,v_p double precision,v_r double precision,q_b double precision);",M);
		else		sprintf(query,"CREATE TABLE master_archive_maps_%s (forecast_time integer,link_id integer,q double precision,s_p double precision,s_l double precision,s_s double precision,v_p double precision,v_r double precision,q_b double precision);",M);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);

		for(i=0;i<numtables;i++)
		{
			if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS archive_maps_%s_%i() INHERITS (master_archive_maps_%s);",M,i,M,M,i);
			else		sprintf(query,"CREATE TABLE archive_maps_%s_%i() INHERITS (master_archive_maps_%s);",M,i,M,M,i);
			res = PQexec(conn,query);
			CheckSQLError(res);
			PQclear(res);

			sprintf(query,"CREATE INDEX idx_archive_maps_%s_%i_forecast_time_link_id ON archive_maps_%s_%i USING btree (forecast_time, link_id);",M,i,M,i);
			res = PQexec(conn,query);
			CheckSQLError(res);
			PQclear(res);
		}
	}

	//Create tables for IFIS

	if(program == 0)
	{
		if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS peakforecast_%s (link_id integer,peak_time integer,peak_discharge double precision,start_time integer,test_stage real);",M);
		else		sprintf(query,"CREATE TABLE peakforecast_%s (link_id integer,peak_time integer,peak_discharge double precision,start_time integer,test_stage real);",M);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);
	}

	if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS hydroforecast_%s(link_id integer NOT NULL,\"time\" integer NOT NULL,discharge double precision,baseflow double precision);",M);
	else		sprintf(query,"CREATE TABLE hydroforecast_%s(link_id integer NOT NULL,\"time\" integer NOT NULL,discharge double precision,baseflow double precision);",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS stageforecast_%s(ifis_id integer,time_utc timestamp with time zone,stage double precision,hr_index integer);",M);
	else		sprintf(query,"CREATE TABLE stageforecast_%s(ifis_id integer,time_utc timestamp with time zone,stage double precision,hr_index integer);",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	if(new_version)	sprintf(query,"CREATE TABLE IF NOT EXISTS warningforecast_%s(ifis_id integer,time_utc timestamp with time zone,warning smallint);",M);
	else		sprintf(query,"CREATE TABLE warningforecast_%s(ifis_id integer,time_utc timestamp with time zone,warning smallint);",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	//Create triggers
	printf("Creating triggers...\n");

	sprintf(query,"CREATE OR REPLACE FUNCTION function_on_insert_to_master_archive_hydroforecast_%s() RETURNS trigger AS $BODY$ BEGIN\
	IF( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= 0        AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < 86400) THEN INSERT INTO archive_hydroforecast_%s_0 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -86400   AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < 0 ) THEN INSERT INTO archive_hydroforecast_%s_1 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -172800  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -86400 ) THEN INSERT INTO archive_hydroforecast_%s_2 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -259200  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -172800) THEN INSERT INTO archive_hydroforecast_%s_3 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -345600  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -259200) THEN INSERT INTO archive_hydroforecast_%s_4 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -432000  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -345600) THEN INSERT INTO archive_hydroforecast_%s_5 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -518400  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -432000) THEN INSERT INTO archive_hydroforecast_%s_6 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -604800  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -518400) THEN INSERT INTO archive_hydroforecast_%s_7 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -691200  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -604800) THEN INSERT INTO archive_hydroforecast_%s_8 VALUES (NEW.*);\
	 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -777600  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -691200) THEN INSERT INTO archive_hydroforecast_%s_9 VALUES (NEW.*);\
	 ELSE RETURN NULL; END IF; RETURN NULL; END; $BODY$\
	LANGUAGE plpgsql VOLATILE COST 100;",M,M,M,M,M,M,M,M,M,M,M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"CREATE TRIGGER trigger_on_insert_to_master_archive_hydroforecast_%s BEFORE INSERT ON master_archive_hydroforecast_%s FOR EACH ROW EXECUTE PROCEDURE function_on_insert_to_master_archive_hydroforecast_%s();",M,M,M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	if(program == 1)
	{
		sprintf(query,"CREATE OR REPLACE FUNCTION function_on_insert_to_master_archive_peakflows_%s() RETURNS trigger AS $BODY$ BEGIN\
		IF( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= 0        AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < 86400) THEN INSERT INTO archive_peakflows_%s_0 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -86400   AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < 0 ) THEN INSERT INTO archive_peakflows_%s_1 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -172800  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -86400 ) THEN INSERT INTO archive_peakflows_%s_2 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -259200  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -172800) THEN INSERT INTO archive_peakflows_%s_3 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -345600  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -259200) THEN INSERT INTO archive_peakflows_%s_4 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -432000  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -345600) THEN INSERT INTO archive_peakflows_%s_5 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -518400  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -432000) THEN INSERT INTO archive_peakflows_%s_6 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -604800  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -518400) THEN INSERT INTO archive_peakflows_%s_7 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -691200  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -604800) THEN INSERT INTO archive_peakflows_%s_8 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -777600  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -691200) THEN INSERT INTO archive_peakflows_%s_9 VALUES (NEW.*);\
		 ELSE RETURN NULL; END IF; RETURN NULL; END; $BODY$\
		LANGUAGE plpgsql VOLATILE COST 100;",M,M,M,M,M,M,M,M,M,M,M);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);

		sprintf(query,"CREATE TRIGGER trigger_on_insert_to_master_archive_peakflows_%s BEFORE INSERT ON master_archive_peakflows_%s FOR EACH ROW EXECUTE PROCEDURE function_on_insert_to_master_archive_peakflows_%s();",M,M,M);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);

		sprintf(query,"CREATE OR REPLACE FUNCTION function_on_insert_to_master_archive_maps_%s() RETURNS trigger AS $BODY$ BEGIN\
		IF( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= 0        AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < 86400) THEN INSERT INTO archive_maps_%s_0 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -86400   AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < 0 ) THEN INSERT INTO archive_maps_%s_1 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -172800  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -86400 ) THEN INSERT INTO archive_maps_%s_2 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -259200  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -172800) THEN INSERT INTO archive_maps_%s_3 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -345600  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -259200) THEN INSERT INTO archive_maps_%s_4 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -432000  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -345600) THEN INSERT INTO archive_maps_%s_5 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -518400  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -432000) THEN INSERT INTO archive_maps_%s_6 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -604800  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -518400) THEN INSERT INTO archive_maps_%s_7 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -691200  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -604800) THEN INSERT INTO archive_maps_%s_8 VALUES (NEW.*);\
		 ELSIF ( new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer >= -777600  AND  new.forecast_time - date_part('epoch', current_date AT TIME ZONE 'UTC')::integer < -691200) THEN INSERT INTO archive_maps_%s_9 VALUES (NEW.*);\
		 ELSE RETURN NULL; END IF; RETURN NULL; END; $BODY$\
		LANGUAGE plpgsql VOLATILE COST 100;",M,M,M,M,M,M,M,M,M,M,M);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);

		sprintf(query,"CREATE TRIGGER trigger_on_insert_to_master_archive_maps_%s BEFORE INSERT ON master_archive_maps_%s FOR EACH ROW EXECUTE PROCEDURE function_on_insert_to_master_archive_maps_%s();",M,M,M);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);
	}

	//Create functions
	printf("Creating functions...\n");

	sprintf(query,"CREATE OR REPLACE FUNCTION copy_to_archive_hydroforecast_%s() RETURNS void AS $BODY$ INSERT INTO master_archive_hydroforecast_%s (link_id,time_utc,discharge,baseflow)\
	(SELECT link_id,to_timestamp(\"time\"),discharge,baseflow FROM hydroforecast_%s);\
	$BODY$ LANGUAGE sql VOLATILE COST 100;",M,M,M);

	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	//For IFIS. These are older functions not used anymore.
/*
	sprintf(query,"CREATE OR REPLACE FUNCTION get_stages_%s()\n\
		  RETURNS void AS\n\
		$BODY$ DROP TABLE IF EXISTS tmp_stageforecast_%s;\n\
		CREATE TABLE tmp_stageforecast_%s AS  (WITH SUPER_I AS (WITH I AS (SELECT ifis_id, max(link_id) as link_id\n\
		   FROM (SELECT * FROM usgs_rating_curves UNION SELECT * FROM ifc_rating_curves) AS K GROUP BY ifis_id)\n\
				SELECT I.*,\n\
				                stage_depth*0.0833333 as last_real,\n\
				                getstages(RR.model_first,I.ifis_id),\n\
				                stage_depth*0.0833333 - getstages(RR.model_first,I.ifis_id) as corr\n\
				FROM I\n\
				LEFT JOIN _link_latest R ON I.ifis_id = R.ifis_id\n\
				LEFT JOIN (SELECT L.link_id, discharge as model_first FROM hydroforecast_%s L INNER JOIN (SELECT link_id, min(\"time\") as min_time FROM hydroforecast_%s group by link_id) as R ON L.link_id = R.link_id AND min_time = \"time\") as RR ON I.link_id = RR.link_id)\n\
		SELECT S_R.ifis_id,\n\
				to_timestamp(S_L.\"time\") as time_utc,\n\
				12*(getstages(S_L.discharge,S_R.ifis_id)+S_R.corr) as stage,\n\
				((\"time\"+1 - EXTRACT('epoch' FROM date_trunc('day', now()) - '10 days'::interval))/3600)::integer as hr_index\n\
		FROM hydroforecast_%s S_L\n\
		INNER JOIN SUPER_I S_R ON S_R.link_id = S_L.link_id);\n\
		\n\
		\n\
		DROP TABLE stageforecast_%s;\n\
		CREATE INDEX ix_stageforecast_%s_ix ON tmp_stageforecast_%s USING btree (ifis_id );\n\
		\n\
		\n\
		ALTER TABLE tmp_stageforecast_%s RENAME TO stageforecast_%s;\n\
		GRANT SELECT ON TABLE stageforecast_%s TO public;\n\
		$BODY$\n\
		  LANGUAGE sql VOLATILE\n\
		  COST 100;"
	,M,M,M,M,M,M,M,M,M,M,M,M,M);

	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"CREATE OR REPLACE FUNCTION update_warnings_%s() RETURNS void AS $BODY$ TRUNCATE warningforecast_%s;\
	INSERT INTO warningforecast_%s (ifis_id,time_utc,warning) (SELECT R.ifis_id,MIN(L.time_utc) AS fcast_t,1 as warning FROM stageforecast_%s L INNER JOIN _link_latest R ON L.ifis_id=R.ifis_id\
	WHERE R.\"case\"=0 AND L.stage>R.action AND R.action>0 and R.action is not null GROUP BY R.ifis_id, R.action);\
	INSERT INTO warningforecast_%s (ifis_id,time_utc,warning) (SELECT R.ifis_id,MIN(L.time_utc) AS fcast_t,1 as warning FROM stageforecast_%s L INNER JOIN _link_latest R ON L.ifis_id=R.ifis_id\
	WHERE R.\"case\"=1 AND L.stage>(R.distance_bottom -R.action) AND R.action>0 and R.action is not null GROUP BY R.ifis_id, R.action);\
	UPDATE warningforecast_%s SET time_utc = fcast_t,warning = 2 FROM (SELECT R.ifis_id,MIN(L.time_utc) AS fcast_t FROM stageforecast_%s L INNER JOIN _link_latest R ON L.ifis_id = R.ifis_id\
	WHERE R.\"case\"=0 AND L.stage > R.flood AND R.flood >0 and R.flood is not null GROUP BY R.ifis_id, R.flood) as I WHERE I.ifis_id=warningforecast_%s.ifis_id;\
	UPDATE warningforecast_%s SET time_utc = fcast_t,warning = 2 FROM (SELECT R.ifis_id,MIN(L.time_utc) AS fcast_t FROM stageforecast_%s L INNER JOIN _link_latest R\
	ON L.ifis_id = R.ifis_id WHERE R.\"case\"=1 AND L.stage > (R.distance_bottom-R.flood) AND R.flood >0 and R.flood is not null GROUP BY R.ifis_id, R.flood) as I WHERE I.ifis_id=warningforecast_%s.ifis_id;\
	UPDATE warningforecast_%s SET time_utc = fcast_t,warning = 3 FROM (SELECT R.ifis_id,MIN(L.time_utc) AS fcast_t FROM stageforecast_%s L INNER JOIN _link_latest R ON L.ifis_id=R.ifis_id\
	WHERE R.\"case\"=0 AND L.stage > R.moderate AND R.moderate>0 and R.moderate is not null GROUP BY R.ifis_id, R.moderate) as I WHERE I.ifis_id=warningforecast_%s.ifis_id;\
	UPDATE warningforecast_%s SET time_utc = fcast_t,warning = 3 FROM (SELECT R.ifis_id,MIN(L.time_utc) AS fcast_t FROM stageforecast_%s L INNER JOIN _link_latest R ON L.ifis_id=R.ifis_id\
	WHERE R.\"case\"=1 AND L.stage > (R.distance_bottom-R.moderate) AND R.moderate>0 and R.moderate is not null GROUP BY R.ifis_id, R.moderate) as I WHERE I.ifis_id=warningforecast_%s.ifis_id;\
	UPDATE warningforecast_%s SET time_utc = fcast_t,warning = 4 FROM(SELECT R.ifis_id, MIN(L.time_utc) AS fcast_t FROM stageforecast_%s L INNER JOIN _link_latest R ON L.ifis_id=R.ifis_id\
	WHERE R.\"case\"=0 AND L.stage > R.major AND R.major>0 and R.major is not null GROUP BY R.ifis_id, R.major) I WHERE I.ifis_id = warningforecast_%s.ifis_id;\
	UPDATE warningforecast_%s SET time_utc = fcast_t,warning = 4 FROM(SELECT R.ifis_id, MIN(L.time_utc) AS fcast_t FROM stageforecast_%s L INNER JOIN _link_latest R ON L.ifis_id=R.ifis_id\
	WHERE R.\"case\"=1 AND L.stage > (R.distance_bottom-R.major) AND R.major>0 and R.major is not null GROUP BY R.ifis_id, R.major) I WHERE I.ifis_id = warningforecast_%s.ifis_id;\
	INSERT INTO warningforecast_%s (ifis_id,time_utc,warning) (SELECT ifisid,time_utc,cast(warninglevel as int) FROM _indexcommunities WHERE cast(warninglevel as int)>0);\
	$BODY$ LANGUAGE sql VOLATILE COST 100;",M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);
*/

	//Clean up
	printf("All done!\n");
	PQfinish(conn);
	return 0;
}

int CheckSQLError(PGresult* res)
{
	short int status = PQresultStatus(res);
	if( !(status == PGRES_COMMAND_OK || status == PGRES_TUPLES_OK) )
	{
		printf("Error making query. %i\n",PQresultStatus(res));
		printf("%s\n",PQresultErrorMessage(res));
		return 1;
	}
	else
		return 0;
}


void CheckConnConnection(PGconn* conn)
{
	int in = PQstatus(conn);
	while(PQstatus(conn) == CONNECTION_BAD)
	{
		printf("Connection to database lost. Attempting to reconnect...\n");
		PQreset(conn);
		sleep(2);
	}
	if(in)	printf("Connection reestablished.\n");
}


//Trashes all tables and functions for a model.
//gcc deletetables.c -o DELETETABLES -O3 -lpq
#include <stdio.h>
#include <stdlib.h>
#include <libpq-fe.h>
#include <string.h>

int CheckSQLError(PGresult* res);

int main(int argc,char* argv[])
{
	int i,j,numtables;
	char query[512];
	char M[32];
	PGresult *res;
	PGconn *conn;

	if(argc < 2)
	{
		printf("Need model name.\n");
		return 1;
	}

	//Connect to database
	//conn = PQconnectdb("dbname=model_ifc host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
	conn = PQconnectdb("dbname=model_test host=s-iihr51.iihr.uiowa.edu port=5432 user=scott password=C5.pfest0");

	//Load up model names
	sprintf(M,argv[1]);
	numtables = 10;
	printf("Deleting tables for model %s.\n\n",M);

	//Delete tables
	for(i=0;i<numtables;i++)
	{
		sprintf(query,"DROP TABLE IF EXISTS archive_hydroforecast_%s_%i;",M,i);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);

		sprintf(query,"DROP TABLE IF EXISTS archive_maps_%s_%i;",M,i);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);

		sprintf(query,"DROP TABLE IF EXISTS archive_peakflows_%s_%i;",M,i);
		res = PQexec(conn,query);
		CheckSQLError(res);
		PQclear(res);
	}

	sprintf(query,"DROP TABLE IF EXISTS master_archive_hydroforecast_%s;",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP TABLE IF EXISTS master_archive_maps_%s;",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP TABLE IF EXISTS master_archive_peakflows_%s;",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP TABLE IF EXISTS hydroforecast_%s;",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP TABLE IF EXISTS peakforecast_%s;",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP TABLE IF EXISTS ratioforecast_%s;",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP TABLE IF EXISTS stageforecast_%s;",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP TABLE IF EXISTS warningforecast_%s;",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP FUNCTION IF EXISTS function_on_insert_to_master_archive_hydroforecast_%s();",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP FUNCTION IF EXISTS function_on_insert_to_master_archive_maps_%s();",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP FUNCTION IF EXISTS function_on_insert_to_master_archive_peakflows_%s();",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP FUNCTION IF EXISTS copy_to_archive_hydroforecast_%s();",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP FUNCTION IF EXISTS get_stages_%s();",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

	sprintf(query,"DROP FUNCTION IF EXISTS update_warnings_%s();",M);
	res = PQexec(conn,query);
	CheckSQLError(res);
	PQclear(res);

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


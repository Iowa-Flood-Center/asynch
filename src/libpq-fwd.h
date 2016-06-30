/* PGconn encapsulates a connection to the backend.
* The contents of this struct are not supposed to be known to applications.
*/
typedef struct pg_conn PGconn;

/* PGresult encapsulates the result of a query (or more precisely, of a single
* SQL command --- a query string given to PQsendQuery can contain multiple
* commands and thus return multiple PGresult objects).
* The contents of this struct are not supposed to be known to applications.
*/
typedef struct pg_result PGresult;

/* PGcancel encapsulates the information needed to cancel a running
* query on an existing connection.
* The contents of this struct are not supposed to be known to applications.
*/
typedef struct pg_cancel PGcancel;
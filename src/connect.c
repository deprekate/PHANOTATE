/* 

    phanotate_connect

    Copyright (C) 2020 Katelyn McNair and Robert Edwards

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


*/

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <getopt.h>
#include <math.h>
//#include "mini-gmp.h"
#include "uthash.h"

#if PY_MAJOR >= 3
#define PY3K
#endif

#define LEN(x) (int)strlen(x)

long double INFINITE = LDBL_MAX;


const char * remove_decimals(const char *str){
	char *output;
	int c, i, len;
	len = 0;
	for (i = 0; i < LEN(str); i++){
		c = (unsigned char) str[i];
		if(c == '.')
			break;
		len = i+1;
	}
	output = malloc(len+1);
	strncpy(output, str, len);
	output[len] = '\0';
	return output;
}

const char * remove_digits(const char *str){
	char* output = malloc(sizeof(char) * LEN(str) );
	int c, i, len;

	len = 0;
	for (i = 0; i < LEN(str); i++){
		c = (unsigned char) str[i];
		if(c < '0' || c > '9'){
			output[len] = c;
			len++;
		}
    }
	output[len] = '\0';
    return output;
}

const char * expand_scinote(const char *str){
	char *output, *ptr;
	size_t int_size, exp_size, size;
	const char *expptr;
	int i, c;

	int_size = 0;
	exp_size = 0;
	//offset = 0;
	expptr = NULL;
	for (i = 0; i < LEN(str); i++){
		c = (unsigned char) str[i];
		//printf("c: %c\n", c);
		if(c == 'e' || c == 'E'){
			if(str[i+1] == '-'){
				output = (char *) malloc(2);
				output[0] = '0';
				output[1] = '\0';
				return output;
			}
			expptr = str + i + 1;
			exp_size = i;
			if(!int_size)
				int_size = i;
			break;
		}else if(c == '.'){
			int_size = i;
		}
	}
	size = int_size + strtol(expptr, &ptr, 10);
	output = malloc(size + 1);
	output[size] = '\0';

	for(i = 0; i < (int)size; i++)
		output[i] = '0';
	for(i = 0; i < (int)int_size; i++)
		output[i] = str[i];
	for(i = int_size+1; i < (int)exp_size; i++)
		output[i-1] = str[i];

	return output;
}

void remove_newline(char *line){
	int new_line = LEN(line) -1;
	if (line[new_line] == '\n')
		line[new_line] = '\0';
}

struct my_struct {
	char key[255];		            /* key */
	char value[255];
	int location;
	char type[255];
	UT_hash_handle hh;		            /* makes this structure hashable */
};

struct my_struct *nodes_left = NULL;	/* important! initialize to NULL */
struct my_struct *nodes_right = NULL;   /* important! initialize to NULL */

void add_leftnode(char *key, char *value) {
	struct my_struct *s;

	HASH_FIND_STR(nodes_left, key, s);  /* id already in the hash? */
	if(s==NULL){
		s = malloc(sizeof(struct my_struct));
		strcpy(s->key, key);
		strcpy(s->value, value);
		s->location = atoi(key);
		strcpy(s->type, remove_digits(key));
		HASH_ADD_STR( nodes_left, key, s );  /* id: name of key field */
	}
}

void add_rightnode(char *key, char *value) {
	struct my_struct *s;

	HASH_FIND_STR(nodes_right, key, s);  /* id already in the hash? */
	if(s==NULL){
		s = malloc(sizeof(struct my_struct));
		strcpy(s->key, key);
		strcpy(s->value, value);
		s->location = atoi(key);
		strcpy(s->type, remove_digits(key));
		HASH_ADD_STR( nodes_right, key, s );  /* id: name of key field */
	}
}

struct not_gaps {
	int location;
	UT_hash_handle hh;		            /* makes this structure hashable */
};
struct not_gaps *iscoding = NULL;	        /* important! initialize to NULL */

void add_coding(int loc) {
	struct not_gaps *s;

	HASH_FIND_INT(iscoding, &loc, s);      /* id already in the hash? */
	if(s==NULL){
		s = malloc(sizeof(struct not_gaps));
		s->location = loc;
		HASH_ADD_INT(iscoding, location, s );  /* id: name of key field */
	}
}

// dont forget to empty

int max_pos = 0;

static PyObject* add_edge (PyObject* self, PyObject* args){
	char *src, *dst, *wgt;
	PyObject *obj;
    if(!PyArg_ParseTuple(args, "O", &obj)) {
        return NULL;
    }
    if(!PyArg_ParseTuple(obj, "sss", &src, &dst, &wgt)) {
        return NULL;
    }
	
	// source
	add_leftnode(src, dst);
	// destination
	add_rightnode(dst, src);

	// do stuff for gap filling
	int i;
	int left = atoi(src);
	int right = atoi(dst);
	for(i=left; i<=right; i++){
		add_coding(i);
	}
	max_pos = (right > max_pos) ? right : max_pos;	

	Py_RETURN_NONE;
}


static PyObject* get_gaps (PyObject* self, PyObject* args, PyObject *kwargs){
	int min_distance = 300;

	static char *kwlist[] = {"min_distance", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|i", kwlist, &min_distance))
		return NULL;

	struct not_gaps *s;
	int *bases = malloc (sizeof (int) * max_pos+1);
	memset (bases, 0, sizeof (int) * max_pos+1);


	for(s=iscoding; s != NULL; s=s->hh.next) {
		bases[s->location] = 1;	
	}

	int i, last = 0, count = 0;
	PyObject *new_edges = PyList_New(0);
	for(i=0; i<=max_pos; i++){
		if(!bases[i]){
			count++;
		}else{
			if(count >= min_distance)
				PyList_Append(new_edges, Py_BuildValue("OOO", PyUnicode_FromFormat("%d_gap", last+1), PyUnicode_FromFormat("%d_gap", i-1), PyUnicode_FromFormat("%d", count) ));
			last = i;
			count = 0;
		}
	}

	return new_edges;
}

static PyObject* get_connections (PyObject* self, PyObject* args, PyObject *kwargs){
	struct my_struct *s1;
    struct my_struct *s2;
	double pnots = 0.99;
	int min_distance = 300;

	static char *kwlist[] = {"pnots", "min_distance", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d|i", kwlist, &pnots, &min_distance))
		return NULL;

	PyObject *new_edges = PyList_New(0);
	int distance;
	// loop over all pairs
	for(s1=nodes_right; s1 != NULL; s1=s1->hh.next) {
			
		for(s2=nodes_left; s2 != NULL; s2=s2->hh.next) {
			// this step is O(n2) so things have to be efficient
			distance = s2->location - (s1->location+3);
			if(-300 < distance && distance < 300){
				distance = (distance >= 0) ? distance/3 : -distance;
				// close by
				if(s1->key != s2->value && s1->value != s2->key){
					//printf("%s - %s\n", s1->key, s2->key);	
					// not the same orf
					if( (atoi(s1->value) < atoi(s2->key)) && (atoi(s2->value) > atoi(s1->key)) ){
						if(strcmp(s1->type, s2->type) !=0 ){
							// same direction
							PyList_Append(new_edges, Py_BuildValue("sss", s1->key, s2->key, PyOS_double_to_string(1/pow(pnots, distance),'f',4,0,NULL) ));
						}else{
							// different direction
							PyList_Append(new_edges, Py_BuildValue("sss", s1->key, s2->key, PyOS_double_to_string(20 + 1/pow(pnots, distance),'f',4,0,NULL) ));
						}
					}
				}
			}
	}
	}
	return new_edges;
}



/*-----------------------------------------------------------------------------------------------*/
/*             Read in node data from stdin                                                      */
/*-----------------------------------------------------------------------------------------------*/

// Our Modules Function Definition struct
// We require this `NULL` to signal the end of our method
static PyMethodDef phanotate_connect_methods[] = {
	{ "get_connections", (PyCFunction) get_connections, METH_VARARGS | METH_KEYWORDS, "Returns the edges of connections orfs" },
	{ "get_gaps",        (PyCFunction) get_gaps,        METH_VARARGS | METH_KEYWORDS, "Returns the gaps as edges" },
	{ "add_edge",        (PyCFunction) add_edge,        METH_VARARGS | METH_KEYWORDS, "Adds an edge to the graph" },
	{ NULL, NULL, 0, NULL }
};
//#ifdef PY3K
// module definition structure for python3
static struct PyModuleDef phanotate_connect = {
	 PyModuleDef_HEAD_INIT,
	"phanotate_connect",
	"mod doc",
	-1,
	phanotate_connect_methods
};
// module initializer for python3
PyMODINIT_FUNC PyInit_connect(void)
{
	return PyModule_Create(&phanotate_connect);
}
/*
#else
 module initializer for python2
PyMODINIT_FUNC initphanotate_connect() {
	Py_InitModule3("phanotate_connect", phanotate_connect_methods, "mod doc");
}
#endif
*/

int main(int argc, char *argv[]){
	wchar_t *program = Py_DecodeLocale(argv[0], NULL);
	if (program == NULL) {
		fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
		exit(1);
	}

	/* Add a built-in module, before Py_Initialize */
	PyImport_AppendInittab("phanotate_connect", PyInit_connect);

	/* Pass argv[0] to the Python interpreter */
	Py_SetProgramName(program);

	/* Initialize the Python interpreter.  Required. */
	Py_Initialize();

	/* Optionally import the module; alternatively,
	   import can be deferred until the embedded script
	   imports it. */
	PyImport_ImportModule("phanotate_connect");

	PyMem_RawFree(program);

	return 0;
}

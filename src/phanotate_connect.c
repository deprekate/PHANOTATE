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
#include "uthash.h"

#if PY_MAJOR >= 3
#define PY3K
#endif

long double INFINITE = LDBL_MAX;


struct my_struct {
    //char key[255];                      /* key */
	int key;
    int value;
    UT_hash_handle hh;                  /* makes this structure hashable */
};

struct my_struct *nodes_left = NULL;    /* important! initialize to NULL */
struct my_struct *nodes_right = NULL;   /* important! initialize to NULL */

void add_leftnode(int key, int value) {
    struct my_struct *s;

    s = malloc(sizeof(struct my_struct));
    //strcpy(s->key, key);
	s->key = key;
    s->value = value;

    HASH_ADD_INT( nodes_left, key, s );  /* id: name of key field */
}

void add_rightnode(int key, int value) {
    struct my_struct *s;

    s = malloc(sizeof(struct my_struct));
    //strcpy(s->key, key);
	s->key = key;
    s->value = value;

    HASH_ADD_INT( nodes_right, key, s );  /* id: name of key field */
}



void remove_newline(char *line){
    int new_line = strlen(line) -1;
    if (line[new_line] == '\n')
        line[new_line] = '\0';
}

static PyObject* add_edge (PyObject* self, PyObject* args, int left_position, int right_position){
	//void add_edge(int edge_id, int src, int dst, long double weight) {
	//char *token;

	//char *edge_string;
	if(!PyArg_ParseTuple(args, "ii", &left_position, &right_position)) {
		return NULL;
	}
	
	// source
	add_leftnode(left_position, right_position);
	// destination
	add_rightnode(right_position, left_position);

	Py_RETURN_NONE;
}

static PyObject* get_connected (PyObject* self, PyObject* args, PyObject *kwargs){
	struct my_struct *s1;
    struct my_struct *s2;
	int min_distance = 300;

	static char *kwlist[] = {"min_distance", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|i", kwlist, &min_distance)) 
	{
		return NULL;
	}

	PyObject *new_edges = PyList_New(0);
	int distance;
    // loop over all pairs
    for(s1=nodes_right; s1 != NULL; s1=s1->hh.next) {
        for(s2=nodes_left; s2 != NULL; s2=s2->hh.next) {
			// this step is O(n2) so things have to be efficient
			distance = s1->key - s2->key;
			distance = (distance >= 0) ? distance : -distance;
			if(distance <= 300 && s1->key != s2->value && s1->value != s2->key){
				PyList_Append(new_edges, Py_BuildValue("iii", s1->key, s2->key, 0));
				//PyList_Append(new_edges, Py_BuildValue("ssi", s1->value, s2->value, 0));
			}
        }
    }
	return new_edges;
}



/*-----------------------------------------------------------------------------------------------*/
/* Read in node data from stdin                                                                  */
/*-----------------------------------------------------------------------------------------------*/
/*
void read_file(){
	struct my_node *node;
	struct my_name *name;
	char buf[256];
	char *token, *err;
	int src, dst;
	long double weight;
	int e = 0;
	int n = 0;
	while (fgets (buf, sizeof(buf), stdin)) {
	}
}
*/

// Our Modules Function Definition struct
// We require this `NULL` to signal the end of our method
static PyMethodDef phanotate_connect_methods[] = {
	{ "get_connected", (PyCFunction) get_connected, METH_VARARGS | METH_KEYWORDS, "Returns the edges of connected orfs" },
	{ "add_edge", (PyCFunction) add_edge, METH_VARARGS | METH_KEYWORDS, "Adds an edge to the graph" },
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
PyMODINIT_FUNC PyInit_phanotate_connect(void)
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

int
main(int argc, char *argv[])
{
	wchar_t *program = Py_DecodeLocale(argv[0], NULL);
	if (program == NULL) {
		fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
		exit(1);
	}

	/* Add a built-in module, before Py_Initialize */
	PyImport_AppendInittab("phanotate_connect", PyInit_phanotate_connect);

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

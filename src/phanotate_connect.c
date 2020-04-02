/* 

    phanotateconnect

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

int n = 0;
int e = 0;

struct my_edge {
	int id;
	int src;
	int dst;
	long double weight;
	UT_hash_handle hh; 
};
struct my_edge *edges = NULL;

struct my_edge *get_edge(int edge_id){
	struct my_edge *s;

	HASH_FIND_INT( edges, &edge_id, s );
	return s;
}

struct my_node {
	char key[256];
	int id;
	UT_hash_handle hh;
};
struct my_node *nodes = NULL;

struct my_name {
	int id;
	char value[256];
	UT_hash_handle hh;
};
struct my_name *names = NULL;


int add_node(char *node_name) {
	struct my_node *node;
	struct my_name *name;
	int src;

	if(!node_name){
		return -1;
	}
	HASH_FIND_STR( nodes, node_name, node);
	if(node == NULL){
		node = (struct my_node*)malloc(sizeof(struct my_node));
		name = (struct my_name*)malloc(sizeof(struct my_name));
		strncpy(node->key, node_name, 256);
		strncpy(name->value, node_name, 256);
		node->id = n;
		name->id = n;
		HASH_ADD_STR( nodes, key, node );
		HASH_ADD_INT( names, id, name );
		src = n; 
		n++;
	}else{
		src = node->id;
	}

	return src;
}

void _add_edge(int edge_id, int src, int dst, long double weight) {
	struct my_edge *s;

	s = (struct my_edge*)malloc(sizeof(struct my_edge));
	s->id = edge_id;
	s->src = src;
	s->dst = dst;
	s->weight = weight;
	HASH_ADD_INT( edges, id, s );
}

static PyObject* add_edge (PyObject* self, PyObject* args){
	//void add_edge(int edge_id, int src, int dst, long double weight) {
	char *token, *err;
	int src, dst;
	long double weight;

	char *edge_string;
	if(!PyArg_ParseTuple(args, "s", &edge_string)) {
		return NULL;
	}

	// parse edge string
	edge_string[strcspn(edge_string, "\n")] = 0;
	token = strtok(edge_string, "\t");	
	// source
	src = add_node(token);
	// destination
	token = strtok(NULL, "\t");
	dst = add_node(token);
	// weight
	token = strtok(NULL, "\t");
	if( src>=0 && dst>=0 && token && (weight = strtold(token, &err)) ){
		_add_edge(e, src, dst, weight);
		e++;
	}else{
		PyErr_SetString(PyExc_ValueError, "Invalid edge");
		return NULL;
	}

	Py_RETURN_NONE;
}

static PyObject* get_connected (PyObject* self, PyObject* args, PyObject *kwargs){
	//(int argc, char *argv[]) {
	struct my_struct *s1;
    struct my_struct *s2;

	char *source, *target, *infile;

	static char *kwlist[] = {(char *)"source", (char *)"target", (char *)"var", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|s", kwlist, &source, &target, &infile)) 
	{
		return NULL;
	}

    // loop over all pairs
    for(s1=nodes_left; s1 != NULL; s1=s1->hh.next) {
        for(s2=nodes_right; s2 != NULL; s2=s2->hh.next) {
            printf("%s\t%s\n", s1->key, s2->key);
        }
    }

	PyObject *path_list = PyList_New(0);
	GetPath(parent, src, dst, path_list);
   
	return path_list;
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
static PyMethodDef phanotateconnect_methods[] = {
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
	phanotateconnect_methods
};
// module initializer for python3
PyMODINIT_FUNC PyInit_phanotateconnect(void)
{
	return PyModule_Create(&phanotate_connect);
}
/*
#else
 module initializer for python2
PyMODINIT_FUNC initphanotateconnect() {
	Py_InitModule3("phanotate_connect", phanotateconnect_methods, "mod doc");
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
	PyImport_AppendInittab("phanotateconnect", PyInit_phanotateconnect);

	/* Pass argv[0] to the Python interpreter */
	Py_SetProgramName(program);

	/* Initialize the Python interpreter.  Required. */
	Py_Initialize();

	/* Optionally import the module; alternatively,
           import can be deferred until the embedded script
           imports it. */
	PyImport_ImportModule("phanotateconnect");

	PyMem_RawFree(program);

	return 0;
}	

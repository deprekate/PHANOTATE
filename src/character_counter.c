static PyObject *
CharCounter(PyObject *self, PyObject *args, PyObject *keywds)
{
    wchar_t *t1;unsigned l1=0;

    if (!PyArg_ParseTuple(args,"u#",&t1,&l1)) return NULL;

    PyObject *resultList,*itemTuple;

    for(unsigned i=0;i<=0xffff;i++)char_counter[i]=0;

    unsigned chlen=0;

    for(unsigned i=0;i<l1;i++){
        if(char_counter[t1[i]]==0)char_list[chlen++]=t1[i];
        char_counter[t1[i]]++;
    }

    resultList = PyList_New(0);

    for(unsigned i=0;i<chlen;i++){
        itemTuple = PyTuple_New(2);

        PyTuple_SetItem(itemTuple, 0,PyUnicode_FromWideChar(&char_list[i],1));
        PyTuple_SetItem(itemTuple, 1,PyInt_FromLong(char_counter[char_list[i]]));

        PyList_Append(resultList, itemTuple);
        Py_DECREF(itemTuple);

    };

    return resultList;
}

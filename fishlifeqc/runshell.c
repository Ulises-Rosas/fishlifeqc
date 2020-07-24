#include <Python.h>
#include <stdio.h>


#define PY_SSIZE_T_CLEAN

// int Cfib(int n){

//     if (n < 2)
//         return n;
//     else
//         return Cfib(n-1) + Cfib(n-2);
// }

static PyObject * get_shell(PyObject *self, PyObject *args){
    const char *command;
    int sts;
    // int n;

    if (!PyArg_ParseTuple(args, "s", &command))
    {
        return NULL;
    }
    // return Py_BuildValue("i", Cfib(n));

    sts = system(command);
    return PyLong_FromLong(sts);
}

static PyObject* version(PyObject* self){

    return Py_BuildValue("s", "Version 0.1");
}

static PyMethodDef runshellmethods[] = {
    {"get", get_shell, METH_VARARGS, "Excute a shell command."},
    {"version", (PyCFunction)version, METH_NOARGS, "Returns the version."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef runshellmodule = {

    PyModuleDef_HEAD_INIT,
    "runshell",
    "Interact with the shell system",
    -1,
    runshellmethods
};

PyMODINIT_FUNC PyInit_runshell(void){

    return PyModule_Create(&runshellmodule);
}

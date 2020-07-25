#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

#define PY_SSIZE_T_CLEAN

// const char* myfile    = "../data/test_hyphen.txt";
// const char* outmyfile = "not_aln_test_hyphen.txt";

// void dealignment(char myarray[], const char* outname);
// char* readfasta(const char* name);

// int main(void){

//     char *myresults = readfasta(myfile);
//     dealignment(myresults, outmyfile); 

//     free(myresults);
// }

char* readfasta(const char* name){

    if (access(name, F_OK) == -1)
    {
        return NULL;
    }

    FILE *f = fopen(name, "rb");
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *string = malloc( fsize + 1 ); // all memory allocation
    if(string){

        fread(string, 1, fsize, f);
        fclose(f);
    }
    else
    {
        fprintf(stderr, "Memory error: Not enough memory to read the file\n");
        fclose(f);
        return NULL;
    }
    string[fsize] = '\0';    
    return string;
}


int dealignment(char myarray[], const char* outname){

    FILE *fp = fopen(outname, "w");
    int itsize = strlen(myarray);

    for (int i = 0; i < itsize; i++)
    {
        // printf("%c", string[i]);
        char mychar = myarray[i];
        if(mychar == '>'){
            do
            {
                fprintf(fp, "%c", mychar);             
                if(mychar == '\n'){
                    break;
                }

                i++;
                mychar = myarray[i];
             
            } while (true);
        }
        else
        {
            if (mychar != '-')
            {
                // mychar = '\0';
                fprintf(fp, "%c", mychar);
            }
        }        
    }

    fprintf(fp, "\n");
    fclose(fp);
    
    return 0;
}

static PyObject* rm_hyphens(PyObject *self, PyObject *args){

    int sts;

    const char* input;
    const char* output;

    if ( !PyArg_ParseTuple(args, "ss", &input, &output) )
    {
        PyErr_SetString(PyExc_SystemError, "Argument parsing issues\n");
        return (PyObject *) NULL;
    }

    char *myresults = readfasta(input);

    if (!myresults)
    {
        PyErr_SetString(PyExc_ValueError, "Issues reading the input file\n");
        return (PyObject *) NULL;
    }

    // dealignment(myresults);
    sts = dealignment(myresults, output);
    free(myresults);
    
    return PyLong_FromLong(sts);
}

static PyObject* version(PyObject* self){
    return Py_BuildValue("s", "version 0.1");
}

static PyMethodDef sequtilsmethods[] = {
    {"rm_hyphens", rm_hyphens, METH_VARARGS, "Remove hyphens from alignment"},
    {"version", (PyCFunction)version, METH_NOARGS, "Returns the version."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef sequtilsmodule = {

    PyModuleDef_HEAD_INIT,
    "fishlifeseq",
    "Sequence utilities",
    -1,
    sequtilsmethods
};

PyMODINIT_FUNC PyInit_fishlifeseq(void){

    return PyModule_Create(&sequtilsmodule);
}


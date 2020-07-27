#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
// #include <unistd.h> // no standar library
#include <math.h>

#define PY_SSIZE_T_CLEAN

// const char* myfile    = "../data/test_hyphen.txt";
// const char* outmyfile = "not_aln_test_hyphen.txt";
// const char* myfile    = "../data/noaln_test_hyphen.txt";

// int dealignment(char myarray[], const char* outname);
// char* readfasta(const char* name);
// int intlen(int number);
// int getheaders(char myarray[], int itsize);
// int get_queries(char myarray[], const char* filename);

// int main(void){

//     char *myresults = readfasta(myfile);
//     get_queries(myresults, myfile);

//     free(myresults);
//     return 0;
// }

char* readfasta(const char* name){

    // if (access(name, F_OK) == -1) // function from <unistd.h>
    // {
    //     return NULL;
    // }

    FILE *f = fopen(name, "rb");

    if (!f)
    {   
        fprintf(stderr, "Couldn't read file\n");
        fclose(f);
        return NULL;
    }

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

int intlen(int number){

    if(!number){
        return 1;
    }

    return floor( log10( abs( number ) ) ) + 1;
}

int getheaders(char myarray[], int itsize){

    // int itsize = strlen(myarray);
    char mychar;
    int headers = 0;

    for (int i = 0; i < itsize; i++)
    {
        mychar = myarray[i];
        if(mychar == '>'){
            headers++;
        }        
    }
    return headers;
}

int get_queries(char myarray[], const char* filename){

    int itsize = strlen(myarray);
    int headers = getheaders(myarray, itsize);

    if (!headers)
    {
        return -1;
    }

    const char* prefix = "_query";
    const int buffer   = strlen(filename) + strlen(prefix);
    char mychar2;

    for (int i = 0; i < headers; i++)
    {
        int   mystrln = buffer + intlen(i) + 1;
        char* OUTFILE = malloc(mystrln);
        snprintf(OUTFILE, mystrln, "%s%s%d", filename, prefix, i);

        FILE *fp      = fopen(OUTFILE, "w");
        int headcount = 0;

        for (int j = 0; j < itsize; j++)
        {
            mychar2 = myarray[j];
            if(mychar2 == '>')
            {
                if (headcount == i)
                {
                    while (true)
                    {
                        fprintf(fp, "%c", mychar2);

                        j++;
                        mychar2 = myarray[j];
                        /*
                        if j is more than itsize,
                        it will overindex and return
                        segmentation fault error
                        */
                        if (mychar2 == '>' || j >= itsize)
                        {
                            break;
                        }
                    }
                    break;
                    /*
                    once matched and extracted all
                    data and break the whole loop 
                    */
                }
                headcount++;
                /*
                count heads until 
                it matchs with `i` with header
                count (sequence number)
                */
            }
        }
        fclose(fp);
    }
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

static PyObject* splitfasta(PyObject *self, PyObject *args){

    int sts;
    const char* input;

    if ( !PyArg_ParseTuple(args, "s", &input) )
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
    sts = get_queries(myresults, input);
    free(myresults);
    
    return PyLong_FromLong(sts);
}

static PyObject* version(PyObject* self){
    return Py_BuildValue("s", "version 0.1");
}

static PyMethodDef sequtilsmethods[] = {
    {"rm_hyphens", rm_hyphens,   METH_VARARGS, "Remove hyphens from alignment"},
    {"splitfasta", splitfasta, METH_VARARGS, "Split fasta by sequence into multiple text files"},
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


#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define PY_SSIZE_T_CLEAN

// const char* myfile    = "../data/test_hyphen.txt";
// const char* outmyfile = "not_aln_test_hyphen.txt";
// const char* myfile    = "../data/noaln_test_hyphen.txt";

// int dealignment(char myarray[], const char* outname);
// char* readfasta(const char* name);
// int intlen(int number);
// int getheaders(char myarray[], int itsize);
// int* mapstrlength(char* myarray, int* tofillup, int itsize, int nheads);
// char** getstrheaders(char* myarray, int nheads);
// int get_queries(char myarray[], const char* filename);

// int main(void){

//     char* filename  = "../data/mock1.txt";
//     char* myresults = readfasta(filename);

//     int itsize = strlen(myresults);
//     int nheads = getheaders(myresults, itsize);

//     char** heads = getstrheaders(myresults, nheads);

//     for (int i = 0; i < nheads; i++)
//     {
//         printf("%s\n", heads[i]);
//     }
    
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

int* mapstrlength(char* myarray, int* tofillup, int itsize, int nheads){

    char mychar;
    int headcount = 0;

    // get length of strings
    for (int i = 0; i < itsize; i++)
    {
        mychar = myarray[i];

        if(mychar == '>'){
            int getterindex = 0;
            do
            { 
                if(mychar == '\n'){
                    break;
                }

                i++;
                getterindex++;

                mychar = myarray[i];
            } while (true);

            tofillup[headcount] = getterindex;
            headcount++;
        } 
    }

    return tofillup;
}


// int** mapstrseqlength(char* myarray, int** tofillup, int itsize, int nheads){

//     char mychar;
//     int headcount = 0;

//     int* headers = tofillup[0];
//     int* seqs    = tofillup[1];

//     // get length of strings
//     for (int i = 0; i < itsize; i++)
//     {
//         mychar = myarray[i];

//         if(mychar == '>'){

//             int headlen = 0;
//             int seqlen  = 0;

//             do
//             { 
//                 if(mychar == '\n'){
//                     break;
//                 }

//                 i++;                
//                 headlen++;
//                 mychar = myarray[i];
//             } while (true);

//             do
//             { 
//                 if(mychar == '>' || i >= itsize){
//                     i--;
//                     break;
//                 }

//                 if(mychar != '\n'){
//                     seqlen++;
//                 }

//                 i++;
//                 mychar = myarray[i];
//             } while (true);

//             tofillup[0][headcount] = headlen;
//             tofillup[1][headcount] = seqlen;
//             headcount++;
//         } 
//     }

//     return tofillup;
// }

char** getstrheaders(char* myarray, int nheads){

    char mychar;

    int itsize = strlen(myarray);
    int strlengths_empty[nheads];

    int*   strlengths = mapstrlength(myarray, strlengths_empty, itsize, nheads);
    char** strheaders = malloc(nheads * sizeof(char*));

    int headcount = 0;
    // get length of strings
    for (int i = 0; i < itsize; i++)
    {
        mychar = myarray[i];

        if(mychar == '>'){
            int getterindex = 0;
            char mytmpstr[ strlengths[headcount] ];
            do
            {                   
                mytmpstr[getterindex] = mychar;

                if(mychar == '\n'){
                    break;
                }

                i++;
                getterindex++;
                mychar = myarray[i];

            } while (true);

            mytmpstr[getterindex] = '\0';

            int mystrln = getterindex + 1;            
            strheaders[headcount] = malloc(mystrln); 
            snprintf(strheaders[headcount], mystrln, "%s", mytmpstr);

            headcount++;
        }
    }
    return strheaders;
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

// static PyObject* concatenate(PyObject *self, PyObject *args){

//     int numLines;
//     char* line;

//     PyObject * listObj; /* the list of strings */
//     PyObject * strObj;  /* one string in the list */

//     if (!PyArg_ParseTuple(args, "O", &listObj))
//     {
//         return NULL;
//     }

//     numLines = PyList_Size(listObj);

//     char** filenames = malloc(numLines * sizeof(char*));
//     int mystrln;

//     for (int i = 0; i < numLines; i++)
//     {
//         strObj         = PyList_GetItem(listObj, i);
//         PyObject* temp = PyUnicode_AsASCIIString(strObj);

//         line = PyBytes_AsString(temp);

//         mystrln      = strlen(line) + 1;
//         filenames[i] = malloc(mystrln);

//         snprintf(filenames[i], mystrln, "%s", line);
//         Py_XDECREF(temp);
//     }

//     for (int i = 0; i < numLines; i++)
//     {
//         printf("%s\n", filenames[i]);
//     }

//     Py_RETURN_NONE;
// }

static PyObject* headers(PyObject *self, PyObject *args){

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

    const int itsize = strlen(myresults);
    const int nheads = getheaders(myresults, itsize);

    if (!nheads)
    {
        PyErr_SetString(PyExc_ValueError, "Check if fasta file is formated properly\n");
        return (PyObject *) NULL;
    }

    char** heads = getstrheaders(myresults, nheads);

    PyObject* python_val = PyList_New(nheads);

    for (int i = 0; i < nheads; ++i)
    {

        PyObject* python_str = Py_BuildValue("s", heads[i]);
        PyList_SetItem(python_val, i, python_str);
        free(heads[i]);
    }

    free(myresults);
    free(heads);

    return python_val;
}

static PyObject* version(PyObject* self){
    return Py_BuildValue("s", "version 0.1");
}

static PyMethodDef sequtilsmethods[] = {
    {"rm_hyphens", rm_hyphens,   METH_VARARGS, "Remove hyphens from alignment"},
    {"splitfasta", splitfasta, METH_VARARGS, "Split fasta by sequence into multiple text files"},
    // {"concatenate", concatenate, METH_VARARGS, "Concatenate fastas"},
    {"headers", headers, METH_VARARGS, "Get sequence headers"},
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


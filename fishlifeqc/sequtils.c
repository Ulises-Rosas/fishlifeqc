#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

const char* myfile = "../E1718.NT_aligned.fasta";

void dealignment(char myarray[]);
char* readfasta(const char* name);

int main(void){

    char *myresults = readfasta(myfile);
    dealignment(myresults);    
    free(myresults);
}

char* readfasta(const char* name){

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
        fprintf(stderr, "Memory error: Not enough memory to read the file");
        fclose(f);
        return NULL;
    }
    string[fsize] = '\0';    
    return string;
}


void dealignment(char myarray[]){

    int itsize = strlen(myarray);

    for (int i = 0; i < itsize; i++)
    {
        // printf("%c", string[i]);
        char mychar = myarray[i];
        if(mychar == '>'){
            do
            {
             printf("%c", mychar);
             
             if(mychar == '\n'){
                 break;
             }

            i++;
            mychar = myarray[i];
             
            } while (true);
        }
        else
        {
            if (mychar == '-')
            {
                mychar = '\0';
            }
            printf("%c", mychar);
        }        
    }
    printf("\n");
}
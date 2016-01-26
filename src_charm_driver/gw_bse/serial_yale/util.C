// Utility functions

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include "util.h"

// Exit program if error occurs.
// simple. no error message. not desired.
void Die(){
    
    printf("\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("                 Program exits.\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    
    exit(EXIT_FAILURE);
}


// print error message
void Die(const char *message){
    
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("%s\n",message);
    printf("Program exits\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    
    exit(EXIT_FAILURE);
}



// print time
void TimeStamp(){
    
    time_t rawtime;       
    struct tm * timeinfo; 

    time(&rawtime); // get current time; same as: timer = time(NULL)
    timeinfo = localtime(&rawtime);

    printf ("Local time: %s", asctime(timeinfo));
    printf (" \n");
}


void mymessage(const char *message){
    printf("\n");
    printf("**  %s\n", message);
}
    

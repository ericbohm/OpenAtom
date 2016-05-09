/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file dictionary.h
 *
 */

#ifndef _Dictionary_
#define _Dictionary_

#define PINY_MAXWORD  80
#define PINY_MAXLINE 100

typedef char PINY_NAME[PINY_MAXWORD];    // Chr: a name; Lth: MAXWORD 
typedef char PINY_LINE[PINY_MAXLINE];    // Chr: a line; Lth: MAXLINE           

typedef struct DICT_WORD{
  int iflag,iuset,key_type;    /* Opt: modifiers                      */
  PINY_NAME keyword,keyarg;    /* Chr: keyword and its arg            */
  PINY_LINE error_mes;         /* Chr: keyword error mess             */
  DICT_WORD(){};
 ~DICT_WORD(){};
}DICT_WORD;

#endif

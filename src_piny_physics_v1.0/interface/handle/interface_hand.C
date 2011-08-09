/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: interface_hand                               */
/*                                                                          */
/* Subprograms that interface with the interface                            */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_handle_entry.h"

/* #define DEBUG */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* dict_print: Function to print dictionary keyword, keyarg, error message  */
/*             based on a type and whether or not the user has set the      */
/*             keyword                                                      */
/*==========================================================================*/

void dict_print(FILE *fp,int num_dict,DICT_WORD dict[],int itype,int iuset)
{ /* begin routine */
  int i;

  for(i=1;i<=num_dict;i++){
    if((dict[i].key_type == itype || itype == -1)  &&
       (dict[i].iuset == iuset || iuset == -1)){
      fprintf(fp,"\\%s{%s}       ERR:%s\n",dict[i].keyword,dict[i].keyarg,
	      dict[i].error_mes);
    }
  }
} /* end routine */

/*==========================================================================*/
/* dict_print_screen: Function to print dictionary keyword (to screen),     */
/*                 keyarg, error message                                    */
/*                 based on a type and whether or not the user has set the  */
/*                 keyword                                                  */
/*==========================================================================*/

void dict_print_screen(int num_dict,DICT_WORD dict[],int itype,int iuset)
{ /* begin routine */
  int i;

  for(i=1;i<=num_dict;i++){
    if((dict[i].key_type == itype || itype == -1)  &&
       (dict[i].iuset == iuset || iuset == -1)){
      PRINTF("\\%s{%s}       ERR:%s\n",dict[i].keyword,dict[i].keyarg,
	      dict[i].error_mes);
    }
  }
} /* end routine */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* get_word: subroutine to read in a keyword and key argument             */
/*==========================================================================*/

int get_word(FILE *fp,DICT_WORD *word,int *nline, int *nkey,
	       int nfun_key,char *file_name)
{
  int ch;
  int nchar_key,nchar_keyarg;
  int backslash,leftcurl,rightcurl,leftbrace,rightbrace,tilde,space,tab;
  int eol;
  int ifind;

  /*========================================================================*/
  /* I) Define identifiers                                                  */

  backslash = (int)'\\';
  leftcurl  = (int)'{';
  rightcurl = (int)'}'; 
  leftbrace = (int)'[';
  rightbrace = (int)']';
  space     = (int)' ';
  tab       = (int)'\t'; 
  tilde     = (int)'~';
  eol       = (int )'\n';

  ifind = 0;  
  nchar_key = 0;
  nchar_keyarg = 0;

  /*========================================================================*/
  /* II) Read until a backslach which is a keyword delimiter or rightbrace */
  
  do{
    ch = fgetc(fp);
    if(ch == eol){(*nline)++;}
    
    if(ch==leftcurl ||ch==rightcurl || ch==leftbrace || ch==tilde || ch==EOF)
      {syntax_error(file_name,*nline,*nkey,nfun_key);}
    
  } while(ch != backslash && ch != rightbrace);
  
  /*========================================================================*/
  /* III) Read keyword (read from backslash to left curly bracket)          */
  
  if(ch == backslash){
    *nkey+=1;
    ifind+=1;
    while( (ch = fgetc(fp)) != leftcurl){

      if(ch == eol){
	(*nline)++;
      }
      if((ch == EOF) || (ch == rightcurl) || (ch == backslash)
	 || (ch == tilde) || (ch == leftbrace) || (ch == rightbrace)){
	PRINTF("ERROR: unexpected character \"%c\"\n",ch);
	syntax_error(file_name,*nline,*nkey,nfun_key);
      }
      if(ch != space && ch != tab){
	if(nchar_key > MAXWORD-2){
	  PRINTF("ERROR: word to long!\n");
	  syntax_error(file_name,*nline,*nkey,nfun_key);
	}
	word->keyword[(nchar_key++)] = (char)ch;
      }
    }/*endwhile*/
    if(nchar_key== 0){
      PRINTF("ERROR:zero length keyword!\n");
      syntax_error(file_name,*nline,*nkey,nfun_key);
    }
    word->keyword[nchar_key] = '\0';
    
    /*=====================================================================*/
    /* C) Read keyarg (read from left curly to right curly bracket)        */
    
    if(ch == leftcurl){
      ifind+=1;
      while( (ch = fgetc(fp)) != rightcurl){
	if(ch == eol){(*nline)++;}
	if((ch == EOF) || (ch == leftcurl) || (ch == backslash)
	   || (ch == tilde) || (ch == leftbrace) || (ch == rightbrace)){
	  PRINTF("ERROR: unexpected character \"%c\"\n",ch);
	  syntax_error(file_name,*nline,*nkey,nfun_key);
	}
	if(ch != space && ch != tab){
	  if(nchar_keyarg+2 == MAXWORD){
	    PRINTF("ERROR: word to long!\n");
	    syntax_error(file_name,*nline,*nkey,nfun_key);
	  }
	  word->keyarg[(nchar_keyarg++)] = (char)ch;
	}
      }
      if(nchar_keyarg== 0){
	PRINTF("ERROR: zero length key arguement\n");
	syntax_error(file_name,*nline,*nkey,nfun_key);
      }
      word->keyarg[nchar_keyarg] = '\0';
    }
  }

  if(ch==rightbrace) return 0;
  else return 1;

} /* end routine */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* syntax_error:  Prints out the syntax error                               */
/*==========================================================================*/

void syntax_error(char *file_name,int nline,int nkey, int nfun_key)
{  /* begin routine */
  PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
  PRINTF("Syntax error occured at or before functional\n");
  if(nfun_key>0){PRINTF("key word number %d\n",nfun_key);}
  if(nkey>0)    {PRINTF("at or before key word number %d\n",nkey);}
  PRINTF("on line %d in file %s\n",nline,file_name);
  PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
  FFLUSH(stdout);
  EXIT(1);
} /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* keyarg_barf:  Pukes when finds an incorrect keyarg                       */
/*==========================================================================*/

void keyarg_barf(DICT_WORD word[],char *file_name, char fun_key[], 
                 int index) 
{ /* begin routine */
  PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
  PRINTF("Bad key argument \"%s\" to key word \\%s\n",
	  word[index].keyarg,word[index].keyword); 
  PRINTF("for functional key word ~%s in file %s\n",
                  fun_key,file_name); 
  PRINTF("Allowed arguments: %s\n",word[index].error_mes);
  PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
  FFLUSH(stdout);
  CkAbort("Bad argument to keyword in input files.");
/* end routine */}
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* keyword_miss:  Pukes when a keyword is not specified                     */
/*==========================================================================*/

void keyword_miss(DICT_WORD word[],char *file_name, char *fun_key,int index)
{ 
  PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
  PRINTF("Required key word \"%s\" to\n",word[index].keyword); 
  PRINTF("functional key word \"%s\" not found in file %s\n",
	  fun_key,file_name); 
  PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
  FFLUSH(stdout);
  EXIT(1);
} 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* put_word_dict: Puts a definition in the dict (overwriting a default def) */
/*==========================================================================*/

void put_word_dict(DICT_WORD *word,DICT_WORD dict[],int num_dict,
                   char fun_key[], int nline,int nkey,int nfun_key,
                   char *file_name)
{ /* begin routine */
  int i,ifound;

  /*=======================================================================*/
  /* I) Puts definition in appropriate location in dictionary              */
  
  ifound = -1;
  i = 0;

  while(++i<=num_dict && ifound == -1){
    if(strcasecmp(word->keyword,dict[i].keyword) == 0){
      strcpy(dict[i].keyarg,word->keyarg);
      dict[i].iuset++;
      ifound = i;
      break;
    }
  }

  /*=======================================================================*/
  /* II) Error checks:                                                     */

  /*----------------------------------------------------------------------*/
  /*  A) Word not found:                                                  */
  
  if(ifound == -1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Bad keyword \"%s\" ",word->keyword);
    if(nfun_key>0){PRINTF("to functional key word \"%s\"\n",fun_key);}
    if(nfun_key>0){PRINTF("at functional key word number %d\n",
                           nfun_key);}
    PRINTF("at key word number %d on line number %d in file %s\n",
	    nkey,nline,file_name);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dictionary and allowed arguments are:\n");
    dict_print_screen(num_dict,dict,-1,-1);
    FFLUSH(stdout);
    EXIT(1);
  }

  /*-----------------------------------------------------------------------*/
  /*  B) Multiple setting of definition:                                   */

  if(dict[ifound].iuset > 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Keyword %s set more then once\n",word->keyword);
    if(nfun_key>0){PRINTF("at functional key word %s\n",fun_key);}
    if(nfun_key>0){PRINTF("at functional key word number %d\n",
                           nfun_key);}
    PRINTF("at key word number %d on line number %d in file %s\n",
	            nkey,nline,file_name);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }
} /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* get_fun_key: Reads in a functional keyword                               */
/*==========================================================================*/
 int get_fun_key(FILE *fp,char fun_key[],int *nline, int *nfun_key,
                 char *file_name)
{ /* begin routine */
/*========================================================================*/
/*             Local variable declarations                                */
/*========================================================================*/
  int ch,nkey;
  int nchar_fun_key;
  int backslash,leftcurl,rightcurl,leftbrace,rightbrace,tilde,space,tab;
  int eol;
/*==========================================================================*/
/* I) Define delimiters:                                                    */

  backslash     = (int)'\\';
  leftcurl      = (int)'{';
  rightcurl     = (int)'}'; 
  leftbrace     = (int)'[';
  rightbrace    = (int)']'; 
  space         = (int)' ';
  tab           = (int)'\t'; 
  tilde         = (int)'~';
  eol           = (int)'\n';
  nchar_fun_key = 0;
  nkey=0;
  
/*==========================================================================*/
/* II) Read until a tilde which is a functional keyword delimiter */
 
  do{
    ch = fgetc(fp);
    if(ch == eol){(*nline)++;}
    
    if(ch==leftcurl ||ch==rightcurl || ch==backslash || ch==leftbrace ||
                      ch==rightbrace)
      {syntax_error(file_name,*nline,nkey,*nfun_key);}
    
  } while(ch != EOF && ch != tilde);

/*==========================================================================*/
/*------------------------------------------------------------------------*/
/* III) Read functional keyword (read from tilde to left brace)             */

  if(ch == tilde){
    *nfun_key+=1;
    while( (ch = fgetc(fp)) != leftbrace){
      if(ch == eol){(*nline)++;}
      if(ch == EOF || ch == rightcurl || ch == leftcurl || ch == backslash
	 || ch == tilde || ch == rightbrace)
	{syntax_error(file_name,*nline,nkey,*nfun_key);}
      if(ch != space && ch != tab)
	{if(nchar_fun_key+2 == MAXWORD)
	   {syntax_error(file_name,*nline,nkey,*nfun_key);}
	 fun_key[(nchar_fun_key++)] = (char)ch;
       }
    }
    if(nchar_fun_key== 0){syntax_error(file_name,*nline,nkey,*nfun_key);}
    fun_key[nchar_fun_key] = '\0';
  }

  if(ch == EOF) return 0;
  else return 1;
} /* end routine */


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* get_fun_key_index: Gets the index of the word in the given dictionary    */
/*                    If word not found, returns -1                         */
/*==========================================================================*/
void get_fun_key_index(char fun_key[],int num_fun_dict,
                      DICT_WORD fun_key_dict[],
                      int nline,int nfun_key,char *file_name,int *num)
{ /* begin routine */
  int i,j;
  
  *num = -1;
  for(i=1;i<=num_fun_dict;i++){
    if(strcasecmp(fun_key_dict[i].keyword,fun_key)==0){
      *num = i;
      fun_key_dict[i].iuset+=1;
      break;
    }
  }
  if(*num==-1){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Functional keyword %s not found \n",fun_key);
     PRINTF("at functional keyword number %d on line %d in file %s\n",
                  nfun_key,nline,file_name);
     PRINTF("Potential keywords are:\n");
     for(j=1;j<=num_fun_dict/2;j++){
       i = j*2-1;
       PRINTF(" %s or ",fun_key_dict[i].keyword);
       PRINTF(" %s or ",fun_key_dict[(i+1)].keyword);
       PRINTF("\n");}
     if((num_fun_dict)%2 !=0){
       PRINTF(" %s ",fun_key_dict[num_fun_dict].keyword);
       PRINTF("\n");
     }
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);}
  i = *num ;
  if(fun_key_dict[i].iuset>1 && fun_key_dict[i].key_type!=2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("ERROR: functional keyword %s set more then once",
                    fun_key_dict[i].keyword);
    PRINTF(" at functional key word number %d on line %d in file %s\n",
	             nfun_key,nline,file_name);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);}

} /* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* get_fun_key_cnt: Reads in a functional keyword                           */
/*==========================================================================*/
int get_fun_key_cnt(FILE *fp,char fun_key[],int *nline, int *nfun_key,
		    char *file_name)
{ /* begin routine */
/*========================================================================*/
/*             Local variable declarations                                */
/*========================================================================*/
  int ch,nkey;
  int nchar_fun_key;
  int backslash,leftcurl,rightcurl,leftbrace,rightbrace,tilde,space,tab;
  int eol;
/*==========================================================================*/
/* I) Define delimiters:                                                    */

  backslash     = (int)'\\';
  leftcurl      = (int)'{';
  rightcurl     = (int)'}'; 
  leftbrace     = (int)'[';
  rightbrace    = (int)']'; 
  space         = (int)' ';
  tab           = (int)'\t'; 
  tilde         = (int)'~';
  eol           = (int)'\n';
  nchar_fun_key = 0;
  nkey=0;
  
/*==========================================================================*/
/* II) Read until a tilde which is a functional keyword delimiter */
 
  do{
    ch = fgetc(fp);
    if(ch == eol){(*nline)++;}
    
    if(ch==leftcurl ||ch==rightcurl || ch==backslash || ch==leftbrace ||
                      ch==rightbrace)
      {syntax_error(file_name,*nline,nkey,*nfun_key);}
    
  } while(ch != EOF && ch != tilde);

/*==========================================================================*/
/* III) Read functional keyword (read from tilde to left brace then         */
/*                             read from left brace to right brace          */
/*                             skipping over key_words and keyargs  )       */
  if(ch == tilde){
    *nfun_key+=1;
    while( (ch = fgetc(fp)) != leftbrace){
      if(ch == eol){(*nline)++;}
      if((ch == EOF) || (ch == rightcurl) || (ch == leftcurl) 
	 || (ch == backslash)
	 || (ch == tilde) || (ch == rightbrace))
	{syntax_error(file_name,*nline,nkey,*nfun_key);}
      if(ch != space && ch != tab)
	{if(nchar_fun_key+2 == MAXWORD)
	   {syntax_error(file_name,*nline,nkey,*nfun_key);}
	 fun_key[(nchar_fun_key++)] = (char)ch;
       }
    }
    if(nchar_fun_key== 0){syntax_error(file_name,*nline,nkey,*nfun_key);}
    fun_key[nchar_fun_key] = '\0';
    while( (ch = fgetc(fp)) != rightbrace){
      if(ch == eol){(*nline)++;}
      if((ch == EOF) || (ch == tilde) || (ch==leftbrace))
	{syntax_error(file_name,*nline,nkey,*nfun_key);}
    }
  } 

  if(ch == EOF) return 0;
  else return 1;
} /* end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* parse_bond_site: Parses keyarg of bond_site into a bond site label       */
/*                  a branch off the bond site label and an                 */
/*                  atom number off the branch label                        */
/*==========================================================================*/
void parse_bond_site(char *strip,char *strip2,NAME site,int *num1,int *num2)
{ 
  int i,j,k;
  strcpy(site,"");
  *num1 = *num2 = -1;
  for(k=0,i=1;i<=3;i++,k++){
    for(j=0;strip[k]!=',' && strip[k]!='\0';j++,k++){
      strip2[j] = strip[k];
    }/*endwhile*/
    strip2[j] = '\0';
    if(i==1&&j>0)sscanf(strip2,"%s",site);
    if(i==2&&j>0)sscanf(strip2,"%d",num1);
    if(i==3&&j>0)sscanf(strip2,"%d",num2);
  }  /*endfor*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* parse_atm_typ_site: Parses keyarg of bond_site into a bond site label    */
/*                  a branch off the bond site label and an                 */
/*                  atom number off the branch label                        */
/*==========================================================================*/
void parse_atm_typ_site(char *strip,char *strip2,int *natm_typ_now,int *iflag,
                        int *count)
{ 
  int i,j,k,iii;
  
/************Fill the vector with the atm types**********************/

  for(i=0;i<1000;i++){
   strip2[i] = '\0';
  }/*endfor*/

  k = *count;j=0;
  while(strip[k]!=',' && strip[k]!='\0'){
   strip2[j] = strip[k];
   j++;k++;
  }/*endwhile*/
   *count = k+1;

   if(strip[k]=='\0')*iflag=1;
   if(strip[k]!=','&&*iflag==1)*iflag=2;

}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* close_fun_key: Reads in a functional keyword                             */
/*==========================================================================*/
void close_fun_key_cnt(FILE *fp,char fun_key[],int *nline, int nfun_key,
		       char *file_name)
{ /* begin routine */
  /*========================================================================*/
  /*             Local variable declarations                                */
  /*========================================================================*/
  int ch,nkey;
  int leftbrace,rightbrace,tilde;
  int eol;
  /*========================================================================*/
  /* I) Define delimiters:                                                  */

  leftbrace     = (int)'[';
  rightbrace    = (int)']'; 
  tilde         = (int)'~';
  eol           = (int)'\n';
  nkey=0;
/*==========================================================================*/
/* II) Reads until functional keyword  closes                               */
  ch = fgetc(fp);
  if(ch != rightbrace){ 
    while( (ch = fgetc(fp)) != rightbrace){
      if(ch == eol){(*nline)++;}
      if((ch == EOF) || (ch == tilde) || (ch == leftbrace))
	{syntax_error(file_name,*nline,nkey,nfun_key);}
    } /*endwhile*/ 
  } /*endif*/
  /*========================================================================*/
}  /* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* dict_save: Function to save dictionary flags and sets and args           */
/*==========================================================================*/
void dict_save(DICT_WORD dict_tmp[],DICT_WORD dict[],int num_dict)
{/* begin routine */
  int i;
  
  for(i=1;i<=num_dict;i++){
    dict[i].iflag    = dict_tmp[i].iflag;
    dict[i].iuset    = dict_tmp[i].iuset;
    sscanf(dict_tmp[i].keyarg,"%s",dict[i].keyarg);    
  }    /*end for*/
}  /* end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* parse_ghost_site: Parses keyarg of bond_site into a bond site label       */
/*                  a branch off the bond site label and an                 */
/*                  atom number off the branch label                        */
/*==========================================================================*/
void parse_ghost(char *strip,char *strip2,int *num1,double *anum2)
{ 
  int i,j,k;
  *num1 = -1;*anum2 = 0.0;
  for(k=0,i=1;i<=2;i++,k++){
    for(j=0;strip[k]!=',' && strip[k]!='\0';j++,k++){
      strip2[j] = strip[k];
    }/*endwhile*/
    strip2[j] = '\0';
    if(i==1&&j>0)sscanf(strip2,"%d",num1);
    if(i==2&&j>0)sscanf(strip2,"%lf",anum2);
  }  /*endfor*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* parse_improp Parses keyarg of improper_ind into 3 integers               */
/*==========================================================================*/
void parse_improp(char *strip,char *strip2,int *num1,int *num2,int *num3,
                                                               int *num4)
{ 
  int i,j,k;
  *num1 = *num2 = *num3 = *num4 = 0;
  for(k=0,i=1;i<=4;i++,k++){
    for(j=0;strip[k]!=',' && strip[k]!='\0';j++,k++){
      strip2[j] = strip[k];
    }/*endwhile*/
    strip2[j] = '\0';
    if(i==1&&j>0)sscanf(strip2,"%d",num1);
    if(i==2&&j>0)sscanf(strip2,"%d",num2);
    if(i==3&&j>0)sscanf(strip2,"%d",num3);
    if(i==4&&j>0)sscanf(strip2,"%d",num4);
  }  /*endfor*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* readtoendofline: Function to read to end of line in read_coord files     */
/*==========================================================================*/
void readtoendofline(FILE *fp){
  int eol,ch;
  eol = (int )'\n';
  ch = eol+1;
  while(ch!=eol){ch=fgetc(fp);}
  if(ch==EOF){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("ERROR: Unexpected end of file reached          \n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
  }/*endif*/  
}/* end routine */    
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* parse_part_lim Parses keyarg of conf_partial_limits into 2 integers      */
/*==========================================================================*/
void parse_part_lim(char *strip,char *strip2,int *num1,int *num2)
{ 
  int i,j,k;
  *num1 = *num2 = 0;
  for(k=0,i=1;i<=2;i++,k++){
    for(j=0;strip[k]!=',' && strip[k]!='\0';j++,k++){
      strip2[j] = strip[k];
    }/*endwhile*/
    strip2[j] = '\0';
    if(i==1&&j>0)sscanf(strip2,"%d",num1);
    if(i==2&&j>0)sscanf(strip2,"%d",num2);
  }  /*endfor*/
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* parse_hydrog_mass Parses keyarg of hydrog_mass_opt with elements  */
/* option and value                                                  */
/*==========================================================================*/
void parse_hydrog_mass(char *strip,char *strip2,int *iopt,double *val,
                       int *ifound)
{
  int i,j,k;
  *iopt = -1;
  *val = 0.0;
  for(k=0,i=1;i<=2;i++,k++){
    for(j=0;strip[k]!=',' && strip[k]!='\0';j++,k++){
      strip2[j] = strip[k];
    }/*endwhile*/
    strip2[j] = '\0';
    if(i==1&&j>0){
      if(strcasecmp(strip2,"off")==0){*iopt=0;}
      if(strcasecmp(strip2,"all")==0){*iopt=1;}
      if(strcasecmp(strip2,"backbone")==0){*iopt=2;}
      if(strcasecmp(strip2,"sidechain")==0){*iopt=3;}
    }
    if(i==2&&j>0)sscanf(strip2,"%lf",val);
  }  /*endfor*/
  *ifound = 1;
  if((*iopt==-1) || (*val<=0.0)){*ifound=0;}
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void parse_on_off(char *keyarg,int *int_val,int *ierr){
  ierr[0]=1;
  if(strcasecmp(keyarg,"off")==0){int_val[0]=0; ierr[0]=0;}
  if(strcasecmp(keyarg,"on")==0) {int_val[0]=1; ierr[0]=0;}
}
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void check_for_slash(char *keyarg,char *keyword, int *ierr_out){
  int len,i,ierr,islash;
  islash = (int )'/';
  len    = (int )(strnlen(keyarg,(size_t)MAXWORD));
  ierr   = 0;
  for(i =0;i<len;i++){
    if(islash==(int ) keyarg[i]){ierr++;}
  }/*end for*/
  if(ierr>0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("No directories permitted in \\%s{%s} : %d detected\n",keyword,keyarg,ierr);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
  }/*endif*/
  ierr_out[0] = ierr;
}
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* readtoendofline_check: Function to read to end of line in input files    */
/* also checks for error, and kicks you out!                                */
/*==========================================================================*/
void readtoendofline_check(FILE *fp,char *f_name,int found, int target){
/*==========================================================================*/

  int ierr,eol,ch;
  eol = (int )'\n';
  ch = eol+1;

/*==========================================================================*/

  while(ch!=eol&&ch!=EOF){ch=fgetc(fp);}
  if(ch==EOF){
    ierr = 1;
  }else{
    ierr = 0;
  }/*endif*/

/*==========================================================================*/

  if(ierr!=0){
    printf("\n@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("something is at odds with this file: %s.\n",f_name);
    if(found<target){
      printf("you said there\'d be %d lines, file ended on line %d.\n",
             target,found);
      printf("Lying will do you no good!!!.\n");
    }/*endif*/
    else{
      printf("methinks you lack or surplus a carriage return\n");
      printf("at the end of line %d\n",found);
    }/*endif*/
    printf("See ya.\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*--------------------------------------------------------------------------*/
 }/* end routine */
/*==========================================================================*/

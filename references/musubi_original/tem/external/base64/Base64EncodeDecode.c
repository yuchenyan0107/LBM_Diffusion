/*
 ============================================================================
 Name        : Base64EncodeDecode.c
 Author      : Sam Ernest Kumar
 Version     :
 Copyright   :
 Description : Base64EncoderDecoder in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Macro definitions */
#define TABLELEN        63

#define ENCODERLEN      4
#define ENCODEROPLEN    0
#define ENCODERBLOCKLEN 3

#define PADDINGCHAR     '='
#define BASE64CHARSET   "ABCDEFGHIJKLMNOPQRSTUVWXYZ"\
                        "abcdefghijklmnopqrstuvwxyz"\
                        "0123456789"\
                        "+/";
/* Function prototypes */
int EncodeIndex(char *input, char *output, int iplen, int index, int ipindex);
int Base64Encode(char *input, char *output, int iplen, int oplen, int ipindex);
int encodeblock(unsigned char *input, unsigned char *output, int iplen, int oplen);
int Base64Decode(char *input, char *output, int oplen);
int decodeblock(char *input, char *output, int oplen);
/* Its always better to move the macros and function prototypes to a header file */

int decodeblock(char *input, char *output, int oplen){
   int rc = 0;
   char decodedstr[ENCODERLEN + 1] = "";

   decodedstr[0] = input[0] << 2 | input[1] >> 4;
   decodedstr[1] = input[1] << 4 | input[2] >> 2;
   decodedstr[2] = input[2] << 6 | input[3] >> 0;
   strncat(output, decodedstr, oplen-strlen(output));

   return rc;
}

int Base64Decode(char *input, char *output, int oplen){
   char *charval = 0;
   char decoderinput[ENCODERLEN + 1] = "";
   char encodingtabe[TABLELEN + 1] = BASE64CHARSET;
   int index = 0, asciival = 0, computeval = 0, iplen = 0, rc = 0;

   iplen = strlen(input);
   while(index < iplen){
      asciival = (int)input[index];
      if(asciival == PADDINGCHAR){
         rc = decodeblock(decoderinput, output, oplen);
         break;
      }else{
         charval = strchr(encodingtabe, asciival);
         if(charval){
            decoderinput[computeval] = charval - encodingtabe;
            computeval = (computeval + 1) % 4;
            if(computeval == 0){
               rc = decodeblock(decoderinput, output, oplen);
               decoderinput[0] = decoderinput[1] =
               decoderinput[2] = decoderinput[3] = 0;
            }
         }
      }
      index++;
   }

   return rc;
}

int encodeblock(unsigned char *input, unsigned char *output, int iplen, int oplen){
   int rc = 0;
   char encodedstr[ENCODERLEN + 2] = "";
   char encodingtabe[TABLELEN + 1] = BASE64CHARSET;

   encodedstr[0] = encodingtabe[ input[0] >> 2 ];
   encodedstr[1] = encodingtabe[ ((input[0] & 0x03) << 4) |
                                 ((input[1] & 0xf0) >> 4) ];
   encodedstr[2] = (iplen > 1 ? encodingtabe[ ((input[1] & 0x0f) << 2) |
                                              ((input[2] & 0xc0) >> 6) ] : PADDINGCHAR);
   encodedstr[3] = (iplen > 2 ? encodingtabe[ input[2] & 0x3f ] : PADDINGCHAR);

   encodedstr[4] = 0;

//   printf("iplen =%d\n", iplen);
//   printf("input 1 = %hhu\n", input[0]);
//   printf("input 2 = %hhu\n", input[1]);
//   printf("input 3 = %hhu\n", input[2]);
//   printf( " encodedstr[0] = %hhu\n", input[0] >> 2 );
//   printf( " encodedstr[1] = %hhu\n", ((input[0] & 0x03) << 4) | ((input[1] & 0xf0) >> 4) );
//   printf( " encodedstr[2] = %hhu\n", ((input[1] & 0x0f) << 2) | ((input[2] & 0xc0) >> 6) );
//   printf( " encodedstr[3] = %hhu\n", input[2] & 0x3f );
//   printf( " encodedstr = %hhu\n", encodedstr);
//   printf( " encodedstr = %hhx\n", encodedstr);

   strncpy(output, encodedstr, ENCODERLEN);

   return rc;
}

int EncodeIndex(char *input, char *output, int iplen, int index, int ipindex){
   int rc = 0;
// printf("input %hhu\n", input[ipindex]);
   if(ipindex < iplen){
      output[index] = input[ipindex];
   }else{
      output[index] = 0;
   }  

   return rc;
}

int Base64Encode(char *input, char *output, int iplen, int oplen, int ipindex){
   int rc = 0;
   int index = 0, ipindex_loc = 0, iplen_sub = 0;
   char encoderinput[ENCODERBLOCKLEN + 1] = "";

   output[0] = 0;
 
   ipindex_loc = ipindex;
   
   for(index = 0; index < 3; index++){
      if(ipindex_loc < iplen){
         encoderinput[index] = input[ipindex];
      }else{
         encoderinput[index] = 0;
      }
      ipindex++;
   }
   iplen_sub = iplen - ipindex_loc + 3;
  
   rc = encodeblock(encoderinput, output, (((iplen_sub) < (3)) ? (iplen_sub) : (3)) , oplen);

//   printf(" output %24s\n", output);
   return rc;
}




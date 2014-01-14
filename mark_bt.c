/* derivative of mark_blast designed to match with .znt files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int,char *[]);

FILE *fp;
FILE *fDatap;

static char chrDtFn[100] = "";
static char chrZntFn[100] = "";
static char chrSet[20] = "";
static char Outp[20] = "btpOut.";
static char znt[5] = ".znt";

main (int argc,char *argv[])
{

        int m, bbChk, readLen, errPrt, bbMax;
        long int cp, chrLen, posMax;
        char chrVal[20], orient, seqRead[120];

        cp = 0;
        chrLen = 0;
        /* readLen = 50; */
        readLen = 40;  /* for tufts_052110 set */

        if(argc != 2)
                {
                printf("argc %d; syntax mark_bt chrNo\n");
                exit(0);
                }

        strcpy(chrSet,argv[1]);

        strcpy(chrDtFn,Outp);
        strcat(chrDtFn,argv[1]);

        strcpy(chrZntFn,chrSet);
        strcat(chrZntFn,znt);

        fp = fopen(chrZntFn,"rb+");
        if (fp == NULL)
                {
                printf("\nOpen chrZntFn rb+ command failed - %s\n",chrZntFn);
                exit(0);
                }
        fseek(fp,0,2);
        chrLen = ftell(fp);
        printf("chrLen  %ld\n",chrLen);

        fDatap = fopen(chrDtFn,"r");
        if (fDatap == NULL)
                {
                printf("\nOpen chrDtFn r command failed - %s\n",chrDtFn);
                exit(0);
                }

        /* mark genome */

        m = 0;
        errPrt = 0;
        bbMax = 0;
        posMax = 0;
        do
                {
                /* fscanf(fDatap," %s %ld %c",chrVal,&cp,&orient); */
                fscanf(fDatap," %s %ld %c %s",chrVal,&cp,&orient,seqRead);  /* seq field added 5/25/2010 */
                if(!strcmp(chrSet,chrVal))
                        {
                        if(m < 5)
                                {
                                printf("chrSet %s, chrVal %s, cp %ld, orient %c\n",chrSet,chrVal,cp,orient);
                                m++;
                                }
                        if(cp >= chrLen)
                                {
                                printf("error - chrLen exceeded\n");
                                exit(0);
                                }

                        if(orient == '+')
                                {
                                fseek(fp,cp - 1,0);
                                bbChk = getc(fp);
                                bbChk = bbChk + 1;

                                if(bbChk <= 128)
                                        {
                                        fseek(fp,cp - 1,0);
                                        putc(bbChk,fp);
                                        }
                                else
                                        {
                                        if(errPrt < 20)
                                                {
                                                printf("error - mark lim 128 exceeded; bbChk %d,  cp %ld\n",bbChk,cp);
                                                errPrt++;
                                                }
                                        }
                                if(bbChk > bbMax)
                                        {
                                        bbMax = bbChk;
                                        posMax = cp;
                                        }
                                }
                        else if(orient == '-')
                                {
                                fseek(fp,cp + readLen - 1,0);
                                bbChk = getc(fp);
                                bbChk = bbChk + 1;

                                if(bbChk <= 128)
                                        {
                                        fseek(fp,cp + readLen - 1,0);
                                        putc(bbChk,fp);
                                        }
                                else
                                        {
                                        if(errPrt < 20)
                                                {
                                                printf("error - mark lim 128 exceeded; bbChk %d,  cp %ld\n",bbChk,cp + readLen);
                                                errPrt++;
                                                }
                                        }
                                if(bbChk > bbMax)
                                        {
                                        bbMax = bbChk;
                                        posMax = cp;
                                        }
                                }
                        }
                }while(!feof(fDatap));

        printf("bbMax %d,  posMax %ld\n\n",bbMax,posMax);

        fclose(fp);
        fclose(fDatap);

        return 0;
}

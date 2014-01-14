/* program to mark annotation on int binary hs genome */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int,char *[]);
int totestMode(void);

FILE *fp;
FILE *fMdatap;
FILE *fTmReadp;
FILE *fTmWritep;

static char chrMdFn[100] = "";
static char chrZntFn[100] = "";
static char tmWFn[100] = "";
static char test[10] = ".test";
static char md[5] = ".md";
static char znt[10] = ".znt";
static char chrChk[10] = "";
static char chr[5] = "chr";

static char gene[10] = "GENE";
static char RNA[5] = "RNA";
static char CDS[5] = "CDS";
static char UTR[5] = "UTR";
static char reference[12] = "reference";
static char Celera[12] = "Celera";

main (int argc,char *argv[])
{

        int nVal, testMode, initChk;
        long int n,chrStart,chrEnd, chrLen;
        char orient, feature[20], group[12];
        char chrSet[20],chrNo[20];

        nVal = 0;
        initChk = 0;
        chrStart = 0;
        chrEnd = 0;
        chrLen = 0;
        orient = '+';

        testMode = 0;

        if(argc != 2)
                {
                printf("argc %d; syntax mark_hsGenome chrNo.znt\n");
                exit(0);
                }

        strcpy(chrSet,argv[1]);

        strcpy(chrMdFn,chrSet);
        strcat(chrMdFn,md);
        strcpy(chrZntFn,chrSet);
        strcat(chrZntFn,znt);
        strcpy(tmWFn,chrSet);
        strcat(tmWFn,test); 

        fp = fopen(chrZntFn,"rb+");
        if (fp == NULL)
                {
                printf("\nOpen chrZntFn rb+ command failed.\n");
                exit(0);
                }
        fseek(fp,0,2);
        chrLen = ftell(fp);
        printf("chrLen  %ld\n",chrLen);

        fMdatap = fopen(chrMdFn,"r");
        if (fMdatap == NULL)
                {
                printf("\nOpen chrMdFn r command failed.\n");
                exit(0);
                }

        /* mark gene regions */

        do
                {
                fscanf(fMdatap," %s %ld %ld %c %s %s",chrNo,&chrStart,&chrEnd,&orient,feature,group);
                if(!initChk)
                        {
                        strcpy(chrChk,chr);
                        strcat(chrChk,chrNo);
                        if(strcmp(chrSet,chrChk))
                                {
                                printf("chrSet/chrChk mismatch\n");
                                printf("chrSet %s; chrChk %s\n",chrSet,chrChk);
                                exit(0);
                                initChk = 1;
                                }
                        }

                if(!strcmp(feature,gene) && orient == '+' && !strcmp(group,reference))
                        {
                        nVal = 1;
                        fseek(fp,chrStart - 1,0);
                        for(n = 0; n < chrEnd - chrStart; n++)
                                {
                                putc(nVal,fp);
                                }
                        }
                if(!strcmp(feature,gene) && orient == '-' && !strcmp(group,reference))
                        {
                        nVal = 2;
                        fseek(fp,chrStart - 1,0);
                        for(n = 0; n < chrEnd - chrStart; n++)
                                {
                                putc(nVal,fp);
                                }
                        }

                }while(!feof(fMdatap));

        /* mark UTR and CDS regions */

        rewind(fMdatap);

        do
                {
                fscanf(fMdatap," %s %ld %ld %c %s %s",chrNo,&chrStart,&chrEnd,&orient,feature,group);

                if((!strcmp(feature,UTR) || !strcmp(feature,CDS)) && orient == '+' && !strcmp(group,reference))
                        {
                        nVal = 3;
                        fseek(fp,chrStart - 1,0);
                        for(n = 0; n < chrEnd - chrStart; n++)
                                {
                                putc(nVal,fp);
                                }
                        }
                if((!strcmp(feature,UTR) || !strcmp(feature,CDS)) && orient == '-' && !strcmp(group,reference))
                        {
                        nVal = 4;
                        fseek(fp,chrStart - 1,0);
                        for(n = 0; n < chrEnd - chrStart; n++)
                                {
                                putc(nVal,fp);
                                }
                        }

                }while(!feof(fMdatap));

        /* mark gene start and end positions */

        rewind(fMdatap);

        do
                {
                fscanf(fMdatap," %s %ld %ld %c %s %s",chrNo,&chrStart,&chrEnd,&orient,feature,group);

                if(!strcmp(feature,gene) && orient == '+' && !strcmp(group,reference))
                        {
                        nVal = 5;
                        fseek(fp,chrStart - 1,0);
                        putc(nVal,fp);

                        nVal = 7;
                        fseek(fp,chrEnd - 1,0);
                        putc(nVal,fp);
                        }
                if(!strcmp(feature,gene) && orient == '-' && !strcmp(group,reference))
                        {
                        nVal = 8;
                        fseek(fp,chrStart - 1,0);
                        putc(nVal,fp);

                        nVal = 6;
                        fseek(fp,chrEnd - 1,0);
                        putc(nVal,fp);
                        }

                }while(!feof(fMdatap));

        fclose(fp);
        fclose(fMdatap);

        if(testMode)
                {
                totestMode();
                }

        printf("done\n");

        return 0;
}

totestMode(void)
{
        int linPtr;
        int bbase;

        fTmReadp = fopen(chrZntFn,"rb");
                if (fTmReadp == NULL)
                {
                printf("chrZntFn rb test command failed.\n");
                exit(0);
                }

        fTmWritep = fopen(tmWFn,"w");
                if (fTmWritep == NULL)
                {
                printf("tmWFn w test command failed.\n");
                exit(0);
                }

        linPtr = 0;

        fprintf(fTmWritep,">%s\n",chrZntFn);
        do
                {
                bbase = getc(fTmReadp);
                fprintf(fTmWritep,"%d",bbase);
                linPtr++;
                if (linPtr == 59)
                        {
                        fprintf(fTmWritep,"\n");
                        linPtr = 0;
                        }
                }while(!feof(fTmReadp));

        fprintf(fTmWritep,"\n");

        fclose(fTmReadp);
        fclose(fTmWritep);

        return 0;
}

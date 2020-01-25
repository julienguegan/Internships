/*
 *  manageData.c
 *
 * This methods are shortcuts to easily manage data.
 *
 * For example, get a value a binary file thanks to the East library (GetInt, GetIntInArray)
 *
 * Or write/read all important data in a file to avoid to have to browse all the binary file and
 * use East library each time of execution (which make loose time because this process can last several minutes)
 *
 *  Created on: Apr 13, 2018
 *  Author: Julien GUEGAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interpreter.h"


/* - Count the number of digit in an integer
 *
 * - input : an int = an integer
 *
 * - output = an int
 */
int LengthInteger(int a_entier) {

	int res = 0;
	if(a_entier == 0) {
		return 1;
	}
	while(a_entier != 0) {
		a_entier /= 10;
		++res;
	}
	return res;
}

/* - Convert a binary number in a relative integer (using the 2's complement)
 *
 * - input : an int = the binary number we want to convert
 *
 * - output : an int
 */
int BinaryToInteger(int a_bin){
	int i = 0, s = 0;

	printf("binaire : %d \n",a_bin);
	int bitDeSigne = a_bin/pow(10,LengthInteger(a_bin)-1);
	printf("bit de signe : %d \n",bitDeSigne);
	if (bitDeSigne == 0) {
		while (a_bin>0){
			s += a_bin%10 * pow(2,i);
			a_bin = a_bin/10;
			i++;}
	} else {
		while (a_bin>0){
			int bit = 0;
			if (a_bin%10 == 0)
				bit = 1;
			s += bit * pow(2,i);
			a_bin = a_bin/10;
			i++;
		}
		s = -s;
	}
	return s;
}

/*  - The error printed is a negative value whose we can check the meaning
 *    by reading the documentation of East library
 *
 *  - input : an int = the output of an East function
 *  - output : if the function didn't work print an error and exit
 */
void CheckEastStatus(int a_status){
	if (a_status == -12) {
		printf("Erreur East : %d mauvais chemin d'accès de la donnée\n", a_status);
		exit(a_status);
	}
	if (a_status != 1) {
		printf("Erreur : %d\n", a_status);
		exit(a_status);
	}
}

void GetLogical(char *a_chaine){
	OCTET *buffer;
	int longueur;
	int resultat = 0;

	resultat = get_data_entity_i(
			a_chaine,
			&buffer,
			&longueur,
			LOGICAL);
	CheckEastStatus(resultat);
}


char* GetAscii(char *a_chaine){
	OCTET *buffer;
	int longueur;
	int resultat = 0;
	char *outputString;

	resultat = get_data_entity_i(
			a_chaine,
			&buffer,
			&longueur,
			ASCII);
	CheckEastStatus(resultat);
	outputString = (char *)buffer;

	free(buffer);

	return outputString;

}

/*  - Used the function get_data_entity() of library East to get an int value, check documentation to know more about it
 *
 *  - input : a string = the acces Path to the asked value define in our model.eas
 *  - output : an int = the value to get back
 */
int GetInt(char *a_chaine){
	OCTET *buffer;
	int longueur;
	int resultat = 0;
	int outputInt = 0;

	resultat = get_data_entity_i(
			a_chaine,
			&buffer,
			&longueur,
			INTEGER_32);
	CheckEastStatus(resultat);
	int *ptr_int = (int *)buffer;
	outputInt = *ptr_int;
	free(buffer);

	return outputInt;

}

/* - Used the function get_data_entity() of library East to get a double value
 *
 */
double GetDouble(char *a_chaine){
	OCTET *buffer;
	int longueur;
	int resultat = 0;
	double outputInt = 0;

	resultat = get_data_entity_i(
			a_chaine,
			&buffer,
			&longueur,
			REAL_64);
	CheckEastStatus(resultat);
	double *ptr_int = (double *)buffer;
	outputInt = *ptr_int;
	free(buffer);

	return outputInt;

}

/* - Concatenates the 2 string to get the full path to data and then call GetInt() function
 *
 * - input : a string = the acces Path (define in our model.eas) to the array which contain a wanted value
 * 			 a string = the remainder of the access Path to go to the specific data contain in the array
 * 			 an int = index of the array toi which we want to recover the data
 * - output : an int = the value to get
 */
int GetIntInArray(char *a_chaineArray, char *a_chaineOutput, int a_index){

	int outputInt = 0;
	char temp[50];
	int taille = strlen(a_chaineArray)+strlen(a_chaineOutput)+LengthInteger(a_index)+3;
	char chaine[taille];

	strcpy(chaine,a_chaineArray);
	sprintf(temp,"(%d).",a_index+1);
	strcat(chaine,temp);
	strcat(chaine,a_chaineOutput);

	outputInt = GetInt(chaine);

	return outputInt;

}

/* - Do the exact same thing of GetIntInArray() but for an array which contains an unique value, we don't need the 2nd string
 *
 * - input : a string = the acces Path (define in our model.eas) to the array which contain a wanted value
 * 			 an int = index of the array toi which we want to recover the data
 * - output : an int = the value to get
 */
int GetIntArray(char *a_chaineArray, int a_index){

	int outputInt = 0;
	char temp[50];
	int taille = strlen(a_chaineArray)+LengthInteger(a_index)+3;
	char chaine[taille];

	strcpy(chaine,a_chaineArray);
	sprintf(temp,"(%d)",a_index+1);
	strcat(chaine,temp);

	outputInt = GetInt(chaine);

	return outputInt;

}

int GetIntIn2Array(char *a_chaineArray1, char *a_chaineArray2, int a_index1, int a_index2, char *a_chaineOutput){

	int outputInt = 0;
	char temp1[5];
	char temp2[5];
	int taille = strlen(a_chaineArray1)+strlen(a_chaineArray2)+strlen(a_chaineOutput)+LengthInteger(a_index1)+LengthInteger(a_index2)+6;
	char chaine[taille];

	strcpy(chaine,a_chaineArray1);
	sprintf(temp1,"(%d).",a_index1+1);
	strcat(chaine,temp1);

	strcat(chaine,a_chaineArray2);
	sprintf(temp2,"(%d).",a_index2+1);
	strcat(chaine,temp2);

	strcat(chaine,a_chaineOutput);

	outputInt = GetInt(chaine);

	return outputInt;

}



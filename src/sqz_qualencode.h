#include <stdint.h>


unsigned char qual_val_table[8] = {33,39,48,55,60,66,70,73};


unsigned char sqz_8binqual(char q);


size_t sqz_qualencode(const char *strqual, uint8_t *codebuff);


size_t sqz_qualdecode(uint8_t *codebuff, char *uncode, size_t length);

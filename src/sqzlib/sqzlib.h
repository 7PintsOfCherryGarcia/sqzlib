/** \file sqzlib.h
Julian Regalado Perez - julian.perez@sund.ku.dk
MIT License

Copyright (c) 2021 Julian Regalado Perez

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>


#define LOAD_SIZE 4L*1024L*1024L   //Sequence buffer size
#define NAME_SIZE 1L*1024L*1024L   //Sequence name buffer size
#define HEADLEN   8UL
#define B64       8UL              //64bits - 8 bytes

/*
  sqzlib file data structure.
  See src/sqz_data.h for details
*/
typedef struct sqzfastx *sqzfastx_t;


/*
 Docstring
*/
typedef struct sqzblock *sqzblock_t;

/*
  Docstring
*/
typedef struct sqzFile_S* sqzFile;

/*
  ##############################################################################
  kseq compatibility routines
  ##############################################################################
*/
/*
  Read len decoded and decompressed bytes into buff
  Returns number of bytes written into buff. If returned value is less than
  len but greater that or equal to 0, end of file has been reached. -1 on error.
*/
int64_t sqzread(sqzFile file, void *buff, uint64_t len);


/*
  Open sqzFile
  Returns sqzFile object or NULL on failure
*/
sqzFile sqzopen(const char *filename, const char *mode);


/*
  Associate an sqzFile with the file descriptor fd
*/
sqzFile sqzdopen(int fd, const char *mode);


/*
 Rewinds sqzFile. Only suported during reading
*/
void sqzrewind(sqzFile file);


/*
  Close sqzFile
*/
void sqzclose(sqzFile file);


/*
##############################################################################
 sqz file routines
##############################################################################
*/

/*
  Load block number n from sqz file
*/
uint8_t sqz_loadblockn(sqzFile sqzfp, uint32_t n);

/*
  Decode up to size bytes of data from sqzFile into buffer
*/
uint64_t sqz_decode(sqzFile sqzfp);

/*
  Write an sqz header
*/
char sqz_filehead(uint8_t fmt, uint8_t libfmt, FILE *ofp);


/*
  Write an sqz tail
*/
uint8_t sqz_filetail(uint64_t numseqs, uint32_t nblocks, FILE *ofp);

/*
  Get format of a fastX file
*/
uint8_t  sqz_getformat(sqzFile sqzfp);


/*
  Get size of an sqz file
*/
uint64_t sqz_filesize(sqzFile sqzfp);

/*
  Get number of blocks in sqz file
*/
uint64_t sqz_getblocks(sqzFile sqzfp);


/*
  Write entire gz sqzFile to ofile
*/
void sqz_gzdump(sqzFile sqzfp, const char *ofile);


/*
  Read len uncompressed bytes from gzip sqzfile into buff
*/
int64_t sqz_gzread(sqzFile sqzfp, void *buff, uint32_t len);


/*
  sqzlib interface to fseek
*/
int sqz_fseek(sqzFile sqzfp, long offset, int whence);


/*
  sqzlib interface to fread
*/
uint64_t sqz_fread(void *ptr, uint64_t size, uint64_t nmemb, sqzFile sqzfp);


/*
  sqzlib interface to ftell
*/
uint64_t sqz_getfilepos(sqzFile sqzfp);


/*
  Returns 1 if sqzfile is fastq format
*/
uint8_t sqz_isfq(sqzFile sqzfp);


/*
  Get compression library
*/
uint8_t sqz_sqzgetcmplib(sqzFile sqzfp);


/*
  Get sqzblock structure from sqzfile
*/
sqzblock_t *sqz_sqzgetblk(sqzFile sqzfp);

/*
  Got to block number in n in sqz file
*/
uint8_t sqz_go2blockn(sqzFile sqzfp, uint64_t n);


/*
  Get format of sqz file
*/
uint8_t sqz_format(sqzFile sqzfp);


/*
  Move decoded block in sqzFile to buffer
*/
uint8_t sqz_emptysqzfp(sqzFile sqzfp, uint8_t *buff);

/*
  ##############################################################################
  sqz structure routines
  ##############################################################################
*/

/*
  Initialize sqzfastx structure
*/
sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t size);


/*
  Free sqzfastx structure
*/
void sqz_fastxkill(sqzfastx_t *sqz);

/*
  Initialize an sqzblock structure
*/
sqzblock_t *sqz_sqzblkinit(uint64_t size);


/*
  Terminate an sqzblock structure
*/
void sqz_blkkill(sqzblock_t *blk);

/*
  Get number of sequences in block
*/
uint64_t sqz_seqsinblk(sqzblock_t *blk);

/*
  ##############################################################################
  sqz data routines
  ##############################################################################
*/

/*
  Read a data block into an sqz block structure
*/
uint8_t sqz_readblksize(sqzFile sqzfp);


/*
  Compress encoded data block
*/
uint64_t sqz_blkcompress(sqzfastx_t *sqz, int level, uint8_t libfmt);


/*
  Encode data from sqz into blk
*/
char sqz_fastXencode(sqzfastx_t *sqz, uint8_t fqflag);


/*
  Decode data in sqz block to memory
*/
uint64_t sqz_fastXdecode(sqzfastx_t *sqz,
                         uint8_t *buff,
                         uint64_t size,
                         uint8_t fqflag);


uint32_t sqz_loadfastX(sqzfastx_t *sqz, sqzFile sqzfp);

/*
  Write compressed block to file
*/
uint8_t sqz_blkdump(sqzfastx_t *sqz, uint64_t cmpsize, FILE *ofp);

/*
  Set block size member to 0
*/
//void sqz_resetblk(sqzblock_t *blk);

/*
  Returns newblk member
*/
uint8_t sqz_newblk(sqzblock_t *blk);

/*
  Decoad and load sqz block into buff
*/
uint64_t sqz_loadblk(sqzblock_t *blk, uint8_t **b, uint64_t *s, uint8_t fqflag);

#ifdef __cplusplus
}
#endif

/*
 * sqzfastx_t data routines
 *
 */

/*
  Check if sqzfastx_t object has data loaded.
  Returns:
      1 if Ture
      0 if False
*/
uint8_t sqz_hasdata(sqzfastx_t *sqz);

/*
  Increase block number and reset appropriate flags
*/
void sqz_resetsqz(sqzfastx_t *sqz);

/*
 Check if reader has ended
*/
uint8_t sqz_readend(sqzfastx_t *sqz);

/*
 Read one last time from sqzfastx_t object
*/
void sqz_setlastread(sqzfastx_t *sqz);

/*
  Set no data flag
*/
void sqz_setnodata(sqzfastx_t *sqz);

/*
  Get number of sequences processed by sqzfastx_t object
*/
uint64_t sqz_getn(sqzfastx_t *sqz);

/** \file sqzlib.h
*/

#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include <stdlib.h>
#include <zlib.h>

#include "sqz_data.h"


/*
    Initialize sqz struct
*/
sqzfastx_t *sqz_fastxinit(const char *filename, uint64_t bsize);


//!SQZLIB file function
/**
 * # Get format of a file
 *
 *
 * Given a file path, sqz_getformat will detect its format. On success
 * a number >= 0 is return corresponding to:
 * 0  - unsupported format
 * 13 - sqz encoded zlib compressed fastA file
 * 14 - sqz encoded zlib compressed fastQ file
 *  2 - fastQ format (compressed or uncompressed)
 *  1 - fastA format (compressed or uncompressed)
 * @param filename
 *
 * @returns unsigned char, Format number
 */
unsigned char sqz_getformat(const char *filename);


/*
    Free sqz object
*/
void sqz_kill(sqzfastx_t *sqz);


/*
    Free blk object
*/
void sqz_killblk(sqzblock_t *blk);
/*
    Initialize blk struct
*/
sqzblock_t *sqz_sqzblkinit(size_t size);

/*
    Write sqz file header
*/
char sqz_filehead(sqzfastx_t *sqz, FILE *ofp);


//!SQZLIB file function
/**
 * # Get size of sqz file excluding header and tail
 *
 *
 * Returns the file size in bytes of an encoded and optionally compressed
 * sqz file. This function fill reset the file position indicator to 0 relative
 * to SEEK_SET (begining of file)
 *
 * @param fp, open file handle
 *
 * @returns uint64_t, Size in bytes of sqzfile
 */
uint64_t sqz_filesize(FILE *fp);


/*
    Load fasta data into sqz struct
*/
size_t sqz_loadfasta(sqzfastx_t *sqz);


/*
    Encode fasta data already loaded in sqz struct and store it in blk struct
*/
char sqz_fastaencode(sqzfastx_t *sqz, sqzblock_t *blk);

/*
    Compress data already loaded in blk struct
*/
size_t sqz_deflate(sqzblock_t *blk, int level);


/*
    Load fastq data into sqz struct
*/
uint64_t sqz_loadfastq(sqzfastx_t *sqz);


/*
    Load fastq data already loaded in sqz struct and store it in blk struct
*/
char sqz_fastqencode(sqzfastx_t *sqz, sqzblock_t *blk);


/*
    Write compressed data stored in blk to file
*/
char sqz_zlibcmpdump(sqzblock_t *blk, size_t size, FILE *ofp);


/*
    Write sqz file tail
*/
char sqz_filetail(size_t numseqs, FILE *ofp);

/*
    Free blk object
*/
void sqz_blkdestroy(sqzblock_t *blk);


/*
    Decode data stored in buff. For fastq data streams
*/
size_t sqz_fastqdecode(sqzblock_t *blk);


/*
    Decode data stored in buff. For fasta data streams
*/
size_t sqz_fastadecode(const uint8_t *buff, size_t size);


//!SQZLIB buffer function
/**
* # Decode sqz data
*
*
* Given blk object with data, an output buffer, and size, sqz_fastXdecode will
* decode up to size decoded bytes and store it in the output buffer.
*
* @param TODO add params
*
* @returns uint64_t, Number of decoded bytes
*/
uint64_t sqz_fastXdecode(sqzblock_t *blk,
                         uint8_t *klibbuff,
                         uint64_t size,
                         char fastqflag);

/*
    Decompress zlib data stored in blk struct
*/
size_t sqz_inflate(sqzblock_t *blk);


/*
#########################################PENE###################################
KSEQ compatibility functions


The reading function must be of the form:
see gzread(): https://github.com/madler/zlib/blob/master/zlib.h
       int reader(FILE *file, void *buff, unsigned len)
Where
    file corresponds to a filepointer
    buff is a non empty buffer where compressed data is extracted
    len number of uncompressed bytes that should be read

This function returns the number of uncompressed bytes read, less than len for
end of file, or -1 for error.

################################################################################
*/



//!kseq compatibility function
/**
 * # Open sqz file for reading
 *
 *
 * Given an sqz filename, sqz_sqzopen opens the file for reading. In this
 * conext, opening a file means initializeg the corresponding structures
 * in an sqz_File object which will be returned by this function. The user
 * is then responsible for closing the sqz_File with sqz_sqzclose().
 *
 * @param filename Path to an sqz file
 *
 * @returns A properly initilized sqz_File object or NULL on error
*/
sqz_File sqz_sqzopen(char *filename);


//! kseq compatibility function
/**
 * # Read uncompressed data in sqz file into a buffer
 *
 *
 * In order to use kseq.h and take advantage of the quick fastX parsing, a
 * "reading" function must be provided. sqz_sqzread loads
 * uncompressed data into a buffer which is subsequently parsed by the
 * rest of kseq.h functions.


 * @param file  sqz_File pointer previously opened with sqz_sqzopen
 * @param buff  Previously allocated buffer of at least "len" bytes
 * @param len   Number of uncomressed bytes to load into "buff"

 * @returns "len" number of loaded bytes, or less than "len" for end of file
            or -1 on error
 */
size_t sqz_sqzread(sqz_File *file, void *buff, size_t len);


//!kseq compatibility function
/**
 * # Close sqz file
 *
 *
 * Given an sqz file, sqz_sqzclose closes a previously opened sqz_File.
 * sqz_sqzclose terminates all necessary objects by calling sqz_kill()
 *
 * @param file sqz_File object previously opened with sqz_sqzopen()
 *
 * @returns Nothing
 */
void sqz_sqzclose(sqz_File file);


//!SQZLIB buffer function
/**
 * # Load and uncompress an sqz block to memory
 *
 * TODO Change all file manipulating functions to deal with sqz_Files instead
        of raw file pointers.
 * This function will read and uncompress an sqz block from file fp.
 * The file position indicator must point to the begining of an sqz block
 * in order to correctly read the sqz block.
 *
 * @param blk, previously initialized sqzblock_t object
 * @param fp,  file pointer
 * @returns 1 on success, 0 on failure
 */
char sqz_readblksize(sqzblock_t *blk, FILE *fp);


//TODO update documentation
//!SQZLIB buffer function
/**
 * # Decodes sqz data block contained in blk which has been previously
 *   decompressed.
 *
 *
 * @param sqz sqzfastx_t object with uncompressed buffer avilable
 * @returns nothing
 */
void sqz_decode(sqzfastx_t *sqz, sqzblock_t *blk, uint64_t klibl);


//!SQZLIB file function
/**
 * # Moves file pointer of sqzfile to the beginning.
 *
 * @param sqzfp - sqzfile to be rewined
 * @return nothing
 */
void sqzrewind(sqz_File *sqzfp);

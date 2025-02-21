
#ifndef _smalloc_h
#define _smalloc_h

/*
 * Memory allocation routines in gromacs:
 *
 * If an allocation fails, the program is halted by means of the
 * fatal_error routine, which outputs source file and line number
 * and the name of the variable involved.
 *
 * Macro's which can be used:
 *
 * snew(ptr,nelem)
 *    Allocates memory for nelem elements and returns this in ptr.
 *    The allocated memory is initialized to zeros.
 *
 * srenew(ptr,nelem)
 *    Reallocates memory for nelem elements and returns this in ptr.
 *
 * smalloc(ptr,size)
 *    Allocates memory for size bytes and returns this in ptr.
 *
 * scalloc(ptr,nelem,elsize)
 *    Allocates memory for nelem elements of size elsize and returns 
 *    this in ptr.
 *
 * srealloc(ptr,size)
 *    Reallocates memory for size bytes and returns this in ptr.
 *
 * sfree(ptr)
 *    Frees memory referenced by ptr.
 *
 ****************************************************************************
 *
 * Functions which are used by the macro's:
 *
 * extern void *save_malloc(char *name,char *file,int line,int size);
 *    Like alloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *
 * extern void *save_calloc(char *name,char *file,int line, 
 *                          unsigned nelem,unsigned elsize);
 *    Like calloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *
 * extern void *save_realloc(char *name,char *file,int line,
 *                           void *ptr,unsigned size);
 *    Like realloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *    If ptr equals NULL, malloc is called in stead of realloc, in this way
 *    it is possible to combine first and later allocations.
 *
 * extern void save_free(char *name,char *file,int line, void *ptr);
 *    Like free, uses name, file and line to generate an error message when 
 *    the free failed.
 *
 * extern unsigned maxavail();
 *    Returns the maximum available allocation unit, by applying a binary
 *    search on the largest block of memory available. After allocation
 *    it invokes free to restore the original state. So it is important
 *    that free can undo the effect of a malloc.
 * 
 * extern unsigned memavail();
 *    Returns the total of available allocation unit, by applying maxavail
 *    until no space is left, it then frees all allocated space and returns
 *    the sum of the previously allocated space. As mentioned with maxavail,
 *    it is important that free can undo the effect of a malloc.
 * 
 */

#define snew(ptr,nelem)   (ptr)=save_calloc((nelem),sizeof(*(ptr)))
#define srenew(ptr,nelem) (ptr)=save_realloc((ptr),(nelem)*sizeof(*(ptr)))
#define smalloc(ptr,size) (ptr)=save_malloc(size)
#define scalloc(ptr,nelem,elsize) (ptr)=save_calloc(nelem,elsize)
#define srealloc(ptr,size)        (ptr)=save_realloc((ptr),size)
#define sfree(ptr)                save_free((ptr))

void *save_malloc(int size); 
void *save_calloc(unsigned nelem,unsigned elsize); 
void *save_realloc(void *ptr,unsigned size);
void save_free(void *ptr);
unsigned maxavail(void);
unsigned memavail(void);

#endif	/* _smalloc_h */

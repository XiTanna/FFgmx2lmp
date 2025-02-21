#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void *save_malloc(int size)
{
  void *p;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
      if ((p=malloc(size))==NULL) 
	{
	  printf("can not allocate memory in save_malloc()\n");
	  exit(0);
	}
      (void) memset(p,0,size);
    }
  return p;
}

void *save_calloc(unsigned nelem,unsigned elsize)
{
  void *p;
  
  p=NULL;
  if ((nelem==0)||(elsize==0))
    p=NULL;
  else
    {
      if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL) 
	{
	  printf("can not allocate memory in save_malloc()\n");
	  exit(0);
	}
    }
  return p;
}

void *save_realloc(void *ptr,unsigned size)
{
  void *p;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
      if (ptr==NULL) 
	p=malloc((size_t)size); 
      else 
	p=realloc(ptr,(size_t)size);
      if (p==NULL) 
	{
	  printf("can not allocate memory in save_malloc()\n");
	  exit(0);
	}
    }
  return p;
}

void save_free(char *name,char *file,int line, void *ptr)
{
  if (ptr != NULL)
    free(ptr);
}

unsigned maxavail(void)
{
  char *ptr;
  unsigned low,high,size;
  
  low=0;
  high=256e6;
  while ((high-low) > 4) {
    size=(high+low)/2;
    if ((ptr=malloc((size_t)size))==NULL)
      high=size;
    else {
      free(ptr);
      low=size;
    }
  }
  return low;
}

unsigned memavail(void)
{
  char *ptr;
  unsigned size;
  
  size = maxavail(); 
  if (size != 0) { 
    if ((ptr=malloc((size_t)size)) != NULL) {
      size += memavail();
      free(ptr);
    }
  }
  return size;
}

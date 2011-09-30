/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: cifplugin.C,v $
 *      $Author:       $       $Locker:  $             $State: Exp $
 *      $Revision:      $       $Date:                     $
 *
 ***************************************************************************/

/* CIF reader plugin */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "molfile_plugin.h"

typedef struct {
  FILE *file;
  int numatoms;
  char *file_name;
  molfile_atom_t *atomlist;
  int coordsloaded;
} cifdata;

  
static void *open_cif_read(const char *filename, const char *filetype, int *natoms) {
  FILE *fd;
  cifdata *data;
  int i;

  fd = fopen(filename, "rb");
  if (!fd) return NULL;

  data = (cifdata *) malloc(sizeof(cifdata));
  memset(data, 0, sizeof(cifdata));

  data->file = fd;
  data->file_name = strdup(filename);
  *natoms = 1;
  data->numatoms=*natoms;
  data->coordsloaded=0;

  return data;
}

static int read_cif_next_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  cifdata *data = (cifdata *)mydata;

  if (data->coordsloaded)
    return MOLFILE_EOF;
	
  for (int i=0; i<natoms; i++) {
    if (ts != NULL) { 
      ts->coords[3*i  ] = 0.012;
      ts->coords[3*i+1] = 0.000;
      ts->coords[3*i+2] = 0.5;
    }
  }

  data->coordsloaded = 1;

  return MOLFILE_SUCCESS;
}


static int read_cif_structure(void *mydata, int *optflags, molfile_atom_t *atoms) {
  *optflags = MOLFILE_NOOPTIONS;

  molfile_atom_t *atom;
  cifdata *data = (cifdata *)mydata;

  for(int i=0; i<data->numatoms; i++) {
    atom = atoms + i;
    strncpy(atom->name, "Se", sizeof(atom->name));
    strncpy(atom->type, "Se1", sizeof(atom->type));
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }

  return MOLFILE_SUCCESS;
}


static int read_cif_bonds(void *v, int *nbonds, int **fromptr, int **toptr, float **bondorderptr) {
  return MOLFILE_SUCCESS;
}

 
static void close_cif_read(void *mydata) {
  cifdata *data = (cifdata *)mydata;
  fclose(data->file);
  free(data->file_name);
  free(data);
}


static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "cif";
  plugin.prettyname = "Crystallographicus";
  plugin.author = "Karel Bartos";
  plugin.majorv = 0;
  plugin.minorv = 1;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "cif";
  plugin.open_file_read = open_cif_read;
  plugin.read_structure = read_cif_structure;
  // plugin.read_bonds = read_cif_bonds;
  plugin.read_next_timestep = read_cif_next_timestep; 
  plugin.close_file_read = close_cif_read;

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}
 
VMDPLUGIN_API int VMDPLUGIN_fini(void) { 
  return VMDPLUGIN_SUCCESS;
}
 



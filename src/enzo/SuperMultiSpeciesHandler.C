/***********************************************************************
/
/  POSTPROCESS CHEMISTRY AND OUTPUT
/
/  written by: Chia-Jung Hsu
/  date:       Sep, 2022j
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include "preincludes.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif
#include "CosmologyParameters.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
                      HierarchyEntry **Grids[]);

int RebuildHierarchy(TopGridData *MetaData,
                     LevelHierarchyEntry *LevelArray[], int level);

int SuperMultiSpeciesHandler(TopGridData *MetaData,
                             LevelHierarchyEntry *LevelArray[], 
                             ExternalBoundary *Exterior
#ifdef TRANSFER
                           , ImplicitProblemABC *ImplicitSolver
#endif
                          )
{

#ifdef USE_NAUNET

  int level;
  LevelHierarchyEntry *Temp;
  HierarchyEntry **Grids;

  /* If we're not using parallel root grid I/O and we're parallel, we
     need to rebuild the hierarchy with the multiple root grids. */

  // if (!ParallelRootGridIO && NumberOfProcessors > 1)
  //   RebuildHierarchy(MetaData, LevelArray, 0);

  /* Compute the potential field on all levels */
  SiblingGridList *SiblingGridListStorage[MAX_DEPTH_OF_HIERARCHY];
  for( int level=0; level < MAX_DEPTH_OF_HIERARCHY; level++ ){
    SiblingGridListStorage[level] = NULL;
  }

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    if (LevelArray[level] != NULL) {
      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
        if (use_naunetstep == 2) {
          Temp->GridData->SyncTopGridCycle(MetaData->CycleNumber);
          Temp->GridData->SetNaunetCycle(MetaData->NaunetCycle);
          Temp->GridData->NaunetWrapper();
          Temp->GridData->SyncNaunetTime();
        }
      } // ENDFOR grids
    } // ENDIF grids on level
  } // ENDFOR levels

#endif

  return SUCCESS;

}

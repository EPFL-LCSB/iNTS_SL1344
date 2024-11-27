# iNTS_SL1344
Genome-scale metabolic model of Salmonella Typhimurium SL1344
and example scripts for the reconstrcuction and analysis performed 
in the manuscript: Vayena et al. "Metabolic network reconstruction as a resource for analyzing
Salmonella Typhimurium SL1344 growth in the mouse intestine"

Requirements
============

You will need to have `Git LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code-block:: shell

git clone https://github.com/EPFL-LCSB/iNTS_SL1344.git
cd iNTS_SL1344
git lfs install
git lfs pull

The scripts have been developed with MATLAB 2019b, and CPLEX 12.7.1 (freely downloadable with the `IBM Academic initiative <https://developer.ibm.com/academic/>`_), and successfully ran on several other versions of both software. However, it is important to respect the IBM compatibility specs sheets between Matlab, CPLEX, and the computer OS - available `on IBM's website <https://www.ibm.com/software/reports/compatibility/clarity/index.html>`_.

This module requires `matTFA <https://github.com/EPFL-LCSB/mattfa>`_, `redGEMX <https://github.com/EPFL-LCSB/redhuman>`_,
`NICEgame <https://github.com/EPFL-LCSB/NICEgame>`_, and `phenomapping <https://github.com/EPFL-LCSB/phenomapping>`_.


For the redGEMX analysis (it requires the redGEMX module):
1. run the script script_prepareModel4redGEMX.m
2. adjust your paths in the functions case_salmonella and get_salmonella in the additionalFunctions folder.
3. run the function get_salmonella


2. Filter the dataset at the restriction fragment level
=======================================================

Before building a heatmap, one have to filter the data at the level of 
individual restriction fragments. This is done via the
:class:`hiclib.fragmentHiC.HiCdataset` class.

In the following code snippet we create a HiCdataset object, populate it with 
the parsed reads generated at the previous step and perform the filtering. 
The following reads are remove from the dataset:

- the reads that start within the 5 bp range from the restriction site
- the identical read pairs, with both ends starting at exactly the same positions
- the reads coming from extremely large and extremely small restriction fragments 
  (length > 10^5 bp or length < 100 bp)
- the reads coming from the top 0.5% most frequently detected restriction fragments

The rationale behind each of the filters is discussed in 
`the hiclib publication <http://dx.doi.org/10.1038/nmeth.2148>`_.
The :doc:`API documentation </fragment>` contains the description of the filters.

.. literalinclude:: ../../../examples/tutorial/02_fragment_filtering.py
   :language: python
   :lines: 1-19

The heatmap is generated from a HiCdataset using the saveHeatmap() function:

.. literalinclude:: ../../../examples/tutorial/02_fragment_filtering.py
   :language: python
   :lines: 20-21

:download:`Download the source code for this chapter. <../../../examples/tutorial/02_fragment_filtering.py>`

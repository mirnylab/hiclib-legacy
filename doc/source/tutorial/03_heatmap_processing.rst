3. Filter and iteratively correct heatmaps.
===========================================

At the last stage of data processing, we filter the data at the binned level 
and perform the iterative correction.

To perform all these operations you have to create 
a :class:`hiclib.binnedData.binnedData` object and load the heatmap that was
generated at the previous step. This class is designed in the way that it can 
store several heatmaps at once, with the same set of filters applied 
to all of them.

.. literalinclude:: ../../../examples/tutorial/03_heatmap_processing.py
   :language: python
   :lines: 1-17

The filters applied at this stage removes the diagonal cells of the heatmap 
as well as the rows and columns corresponding to the poorly sequenced or covered
regions of the genome, and the tentative PCR blowouts. For the detailed
discussion of the filters and the order that they should be applied, we refer
the reader to :mod:`the API documentation <hiclib.binnedData>`.

.. literalinclude:: ../../../examples/tutorial/03_heatmap_processing.py
   :language: python
   :lines: 18-32

After all the corrections, you can either access the heatmap directly through
the dataDict member variable of :class:`hiclib.binnedData.binnedData` or
you can export it into an h5dict file.

.. literalinclude:: ../../../examples/tutorial/03_heatmap_processing.py
   :language: python
   :lines: 33-38

:download:`Download the source code for this chapter. <../../../examples/tutorial/03_heatmap_processing.py>`

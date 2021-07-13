def plot_neighboreffect(self, offset=1, output_file: Union[None, str, Path] = None, **kwargs):
    """
   DEPRECATED.

   Generate a miniheatmap plot telling you the effect of having a
   residue in front or behind. It corrects for the effect of that amino
   acid on the rest of the population.

   Parameters
   ----------
   self : object from class *Screen*

   offset : int, default 1
      if you want to study effects of a residue when is behind or in
      front of another residue. offset of 1 means that you evaluate the
      effect of following residue n+1 on n. On a "MTEY..." sequence,
      you would look at the effect of T on M, E on T, Y on E, etc.. and
      then group by residue (n+1). offset of -1 means that you look at
      the previous residue (n-1 on n).

   output_file : str, default None
      If you want to export the generated graph, add the path and name
      of the file. Example: 'path/filename.png' or 'path/filename.svg'.

   **kwargs : other keyword arguments
       return_plot_object : boolean, default False
           If true, will return plotting objects (ie. fig, ax).

   Returns
   ----------
   fig, ax, cb : matplotlib figure and subplots
       Needs to have return_plot_object==True. By default they do
       not get returned.

   """
    # load font parameters
    code_kwargs._font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    if '*' in temp_kwargs['neworder_aminoacids']:
        temp_kwargs['neworder_aminoacids'].remove('*')

    # do offset, no stop codons
    df = _normalize_neighboreffect(self, offset, temp_kwargs['neworder_aminoacids'])

    # Plot
    fig, ax, cb = _plot_miniheatmap(df, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax, cb

    # show figure
    if temp_kwargs['show']:
        plt.show()

pro plot_rectangular_som, somfile

readcol, somfile, index, X, Y, Weight1, Weight2, Weight3, Weight4, Weight5, Weight6, Weight7, Weight8, Cell_Occupation, U_value, Cell_Dispersion, Cell_Average_Redshift
 
; Set up variables for the plot. Normally, these values would be 
; passed into the program as positional and keyword parameters.
;x = cgScaleVector(Randomn(-3L, 100000)*2, -10, 10)
;y = cgScaleVector(Randomn(-5L, 100000)*5, 0, 100)
;xrange = [Min(x), Max(x)]
;yrange = [Min(y), Max(y)]
;xbinsize = 0.25
;ybinsize = 3.00

xrange = [min(X), max(X)]
yrange = [min(Y), max(Y)]

cell_occupation_array = dblarr(max(X)+1,max(Y)+1)
cell_u_array = dblarr(max(X)+1,max(Y)+1)
cell_disp_array = dblarr(max(X)+1,max(Y)+1)
cell_redshift_array = dblarr(max(X)+1,max(Y)+1)

for i = 0, max(X) do begin
   for j = 0, max(Y) do begin
      thisindex = j*(max(X)+1) + i; check this
      cell_occupation_array[i,j] = Cell_Occupation[thisindex]
      cell_u_array[i,j] = U_value[thisindex]
      cell_disp_array[i,j] = Cell_Dispersion[thisindex]
      cell_redshift_array[i,j] = Cell_Average_Redshift[thisindex]
   endfor
endfor
      

; Open a display window.
cgDisplay

; Create the density plot by binning the data into a 2D histogram.
;density = Hist_2D(x, y, Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize, $
;                        Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize)   
                        
;maxOccupation = Ceil(Max(cell_occupation_array)/1e2) * 1e2
scaledOccupationArray = BytScl(cell_occupation_array, Min=0, Max=max(cell_occupation_array))

;maxU = Ceil(Max(cell_u_array)/1e2) * 1e2
scaledUArray = BytScl(cell_u_array, Min=0, Max=max(cell_u_array))

;maxDisp = Ceil(Max(cell_disp_array)/1e2) * 1e2
scaledDispArray = BytScl(cell_disp_array, Min=0, Max=max(cell_disp_array))

;maxRedshift = Ceil(Max(cell_redshift_array)/1e2) * 1e2
scaledRedshiftArray = BytScl(cell_redshift_array, Min=0, Max=max(cell_redshift_array))
                      
; Load the color table for the display. All zero values will be gray.
cgLoadCT, 33
TVLCT, cgColor('gray', /Triple), 0
TVLCT, r, g, b, /Get
palette = [ [r], [g], [b] ]

 cgDisplay, 1200, 1200
   pos = cgLayout([2,2], OXMargin=[3, 3], OYMargin=[3, 14], XGap=5, YGap=12)

   ;cgLoadCT, 3, /Reverse, /Brewer, RGB_Table=palette
   ;cgImage, cgDemoData(18), Position=pos[*,0], Palette=palette, /Save
   cgImage, scaledRedshiftArray, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
   position=pos[*,0], /Save
   cbpos = [!X.Window[0], !Y.Window[1]+0.05, !X.Window[1], !Y.Window[1]+0.10]
   ;cgColorbar, Palette=palette, Position=cbpos
   cgColorbar, Position=cbpos, color='000000'XL, Title='Mean Redshift', Palette=palette, $
    Range=[0, max(cell_redshift_array)], NColors=254, Bottom=1, OOB_Low='gray'


   cgLoadCT, 5, /Reverse, /Brewer, RGB_Table=palette
   ;cgImage, cgDemoData(18), Palette=palette, /Save, /NoErase
   cgImage, scaledUArray, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
   Position=pos[*,1], /Save, /noerase
   cbpos = [!X.Window[0], !Y.Window[1]+0.05, !X.Window[1], !Y.Window[1]+0.10]
   ;cgColorbar, Palette=palette, Position=cbpos
   cgColorbar, Position=cbpos, color='000000'XL, Title='U Value', Palette=palette, $
    Range=[0, max(cell_u_array)], NColors=254, Bottom=1, OOB_Low='gray'

   cgLoadCT, 7, /Brewer, RGB_Table=palette
   ;cgImage, cgDemoData(18), Palette=palette, /Save, /NoErase
   cgImage, scaledDispArray, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
   Position=pos[*,2], /Save, /noerase
   cbpos = [!X.Window[0], !Y.Window[1]+0.05, !X.Window[1], !Y.Window[1]+0.10]
   ;cgColorbar, Palette=palette, Position=cbpos
   cgColorbar, Position=cbpos, color='000000'XL, Title='Cell Dispersion', Palette=palette, $
    Range=[0, max(cell_disp_array)], NColors=254, Bottom=1, OOB_Low='gray'

  ; cgLoadCT, 7, /Reverse, /Brewer, RGB_Table=palette
  ; cgImage, cgDemoData(18), Position=pos[*,2], Palette=palette, /Save, /NoErase
  ; cbpos = [!X.Window[0], !Y.Window[1]+0.05, !X.Window[1], !Y.Window[1]+0.10]
  ; cgColorbar, Palette=palette, Position=cbpos

   cgLoadCT, 9, /Brewer, RGB_Table=palette
   ;cgImage, cgDemoData(18), Palette=palette, /Save, /NoErase
   cgImage, scaledOccupationArray, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
   Position=pos[*,3], /Save, /noerase
   cbpos = [!X.Window[0], !Y.Window[1]+0.05, !X.Window[1], !Y.Window[1]+0.10]
   ;cgColorbar, Palette=palette, Position=cbpos
   cgColorbar, Position=cbpos, color='000000'XL, Title='Cell Occupation', Palette=palette, $
    Range=[0, max(cell_occupation_array)], NColors=254, Bottom=1, OOB_Low='gray'

 cgText, 0.5, 0.97, Alignment=0.5, /Normal, 'SOM Overview', Charsize=3, Charthick=2, color='000000'XL

  ; cgLoadCT, 9, /Reverse, /Brewer, RGB_Table=palette
  ; cgImage, cgDemoData(18), Position=pos[*,3], Palette=palette, /Save, /NoErase
  ; cbpos = [!X.Window[0], !Y.Window[1]+0.05, !X.Window[1], !Y.Window[1]+0.10]
  ; cgColorbar, Palette=palette, Position=cbpos
  ; cgText, 0.5, 0.95, Alignment=0.5, /Normal, 'This is the Plot Title', Charsize=1.75

end ;*****************************************************************

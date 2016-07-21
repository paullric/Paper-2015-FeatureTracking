;***********************************************
; plot-etc-density.ncl
;***********************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

;************************************************
begin

;************************************************
; read in netCDF file
;************************************************
  density_file = addfile("cfsr/pgbhnl.dcu_etc_density.nc","r")

  density = density_file->density

;************************************************
; plot resources [options]
;************************************************
  plot = new(1,graphic)                ; define a panel plot

; Resources:
  res                     = True

  res@gsnDraw             = False
  res@gsnFrame            = False

  res@vpKeepAspect        = True
  res@tiMainFontHeightF   = 0.02
  res@vpWidthF            = 0.8
  res@vpHeightF           = 0.4
;  res@gsnMaximize         = True
;  res@gsnPaperOrientation = "portrait"

  wks = gsn_open_wks("eps","plot-etc-density")   

; choose color table:
;  gsn_define_colormap(wks,"cosam12")             ; not enough colors
;  gsn_define_colormap(wks,"WhBlGrYeRe")
;  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")   ; favorite right now
;  gsn_define_colormap(wks,"BlGrYeOrReVi200")
;  gsn_define_colormap(wks,"BlueRed")              ; mid-light blue to red
;  gsn_define_colormap(wks,"BlRe")                ; too dark
;  gsn_define_colormap(wks,"wh-bl-gr-ye-re")
;  gsn_define_colormap(wks,"WhBlGrYeRe")
;  gsn_merge_colormaps(wks,"gui_default","BlueDarkRed18")
  gsn_define_colormap(wks,"gui_default")

  colors = gsn_retrieve_colormap(wks)
  colors(2,0) = 1.0
  colors(2,1) = 1.0
  colors(2,2) = 1.0
  gsn_define_colormap(wks, colors)

  res@cnFillOn          = True
  res@gsnSpreadColors   = True

  res@cnLinesOn         = False
  res@cnLineLabelsOn    = False
  res@lbLabelAutoStride = True
  res@cnInfoLabelOn     = False
  res@tmXBLabelsOn      = True
  res@tmXBOn            = True

  res@tmXTLabelsOn      = False
  res@tmXTOn            = False

  res@gsnContourZeroLineThicknessF = 0
  res@gsnContourNegLineDashPattern = 6

  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200/)

  res@lbLabelBarOn = True
  res@lbLabelFontHeightF = 0.018

  res@tiMainFontHeightF = 0.018

  res@gsnSpreadColorStart = 2
  res@gsnSpreadColorEnd = 23

  res@tiMainString = "Cold-Core Cyclone Counts from CFSR"
  plot(0) = gsn_csm_contour_map(wks, density, res)

  pres = True
  gsn_panel(wks, plot, (/1,1/), pres)

end

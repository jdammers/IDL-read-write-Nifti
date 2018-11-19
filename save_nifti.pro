;; ====================================================================================
;;
;;  save_nifti
;;
;; ARGUMENTS
;;   fnifti:    output filename of nifti header
;;   vol   :    3D volume to save (optional)
;;
;; OPTIONS
;;   header:  use this header as a header template 
;;            (otherwise defaults are used as stored in IDLffNifti__define.pro)
;;
;;   just_header: store only header file, no volume will be stored 
;;       - dimension: dimension of the volume (used togehter with just_header)
;;       - byte_type: idl byte type of the data (used togehter with just_header)
;;
;;   show:    show header info on screen after saving the data
;;
;;   dimension:  user defined image/volume dimension
;;   pixdim:  dimension of the voxel in mm
;;            default is pixdim = [1.0,1.0,1.0]
;;
;;   origin:  set new origin 
;;
;;   vox_offset:  set new voxel offset
;;   quatern_bcd: set quaterions b,c,d
;;   smatrix    : instead of quatern_bcd you can provide the full 3x4 srow matrix
;;                example:
;;                smatrix = [ $
;;                          [ [-1, 0, 0] * pixdim[0], origin[0] ], $ ; left - right  orientation (radiological)
;;                          [ [ 0, 0,-1] * pixdim[1], origin[1] ], $ ; inferior - superior orientation  (from bottom to top)
;;                          [ [ 0, 1, 0] * pixdim[2], origin[2] ]]   ; posterior - anterior
;;
;;   noreverse:  do not flip top-bottom image orientation
;;               ==> the default is to flip the image because of IDL's image origin
;;
;; requirements:
;;    - needs:  IDLffNifti__define.pro
;;
;; Author: Juergen Dammers, J.Dammers@fz-juelich.de
;; ====================================================================================
FUNCTION save_nifti, fnifti, vol, header=hdr, just_header=just_header, dimension=dimension, show=show, $
         pixdim=pixdim, origin=origin, vox_offset=vox_offset,quatern_bcd=quatern_bcd, $
         noreverse=noreverse, smatrix=smatrix


  ;; ------------------------------------------------------------------
  ;; check input
  ;; ------------------------------------------------------------------
  sz  = size(vol)
  sz_dim = size(vol,/dim)
  if (sz_dim[0] ne 0) then ndim = n_elements(sz_dim) else ndim=0
  if (total(sz) lt 3) then just_header=1 
  ext   = strmid(fnifti,strlen(fnifti)-3,3)
  if (ext eq 'nii') then begin
      if not keyword_set(vox_offset) then vox_offset=352
  endif else vox_offset=0
  if (ndim eq 4) then nvol = sz_dim[3] else nvol=1
  

  ;; ------------------------------------------------------------------
  ;; create Nifti Object
  ;; ------------------------------------------------------------------
  oNifti = obj_new('IDLffNifti') 



  ;; ------------------------------------------------------------------
  ;; pass input header
  ;; ------------------------------------------------------------------
  if keyword_set(hdr) then begin
      header = hdr
      oNifti -> SetProperty,hdr 
  endif else header = oNifti -> GetProperty()


  ;; ------------------------------------------------------------------
  ;; check keywords
  ;; ------------------------------------------------------------------    
  ;; dimension
  if keyword_set(dimension) then dim=dimension else begin
     if (ndim gt 1) then begin
        dim         = make_array(8,/integer,value=1)
        dim[0]      = ndim
        dim[1:ndim] = sz_dim[0:ndim-1]
        dim[4]      = nvol
     endif else dim = header.dim
  endelse
  dim_ok = onifti->isnum(dim,/Positive) * ndim
  if (dim_ok lt ndim) then begin
     tmp = dialog_message('ERROR (save_nifti): Invalid dimension: must be 3, 4 or 5.')
     return,-1
  endif 
  oNifti -> SetDim, dim

  ;; origin
  if keyword_set(origin) then begin
     orig_ok = (onifti->isnum(origin)) * n_elements(origin)
     origin  =  float(origin)
     if (orig_ok eq 3) then oNifti -> SetOrigin, origin
  endif

  ;; pixdim
  if keyword_set(pixdim) then begin
     ok = (onifti->isnum(pixdim)) * n_elements(pixdim)
     pixdim  =  float(pixdim)
     pixdim_array = header.pixdim
     pixdim_array[1:3] = pixdim
     if (ok eq 3) then oNifti -> SetPixdim, pixdim_array
  endif

  ;; sform matrix
  if keyword_set(smatrix) then begin
     oNifti->CalcQuatern,smatrix=smatrix
  endif
  
   ;; quatern_bcd
  if keyword_set(quatern_bcd) then oNifti->SetQuatern_bcd, float(quatern_bcd)

  ;; extension
  if (ext eq 'nii') and (header.vox_offset eq 0) then oNifti -> SetVoxoffset,vox_offset



  ;; ------------------------------------------------------------------
  ;; pass data
  ;; ------------------------------------------------------------------    
  if not keyword_set(just_header) then begin
     ;; check if we need to flip top -> bottom
      if not keyword_set(noreverse) then volume = reverse(vol,2) else volume=vol
      oNifti -> SetData, volume
      header =  oNifti -> GetProperty()
  endif



  ;; ------------------------------------------------------------------
  ;; write data
  ;; and get header on return
  ;; ------------------------------------------------------------------
  oNifti -> WriteFile, fnifti, just_header=just_header
  header = oNifti -> GetProperty()


  ;; ------------------------------------------------------------------
  ;; show header info
  ;; ------------------------------------------------------------------
  if keyword_set(show) then begin
     print, '>>> filename:'
     print, fnifti
     oNifti->OutputHeader, header,-1
  endif
  
  obj_destroy, oNifti

  return,header
END

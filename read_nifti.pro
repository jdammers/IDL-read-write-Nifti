;; ====================================================================================
;; read_nifti
;;
;; ARGUMENTS
;; fname        Nifti filename (either .nii or a pair of .hdr and .img)
;;
;;
;; OPTIONS      
;; header      to get a structure with the Nifti header on return
;;
;; just_header reads the header information only
;; 
;; show        show header information on screen
;;
;; compress    reads gnuzip Nifti files
;;
;; flip_ver    performs a top/bottom flip on the image/volume
;;             ==> Note this is the default
;;
;; flip_hor    performs a horizontal flip on the volume (e.g. for HyperView)
;;             
;;
;; requirements:
;;    - needs:  IDLffNifti__define.pro
;; 
;; Author: Juergen Dammers, J.Dammers@fz-juelich.de
;; ====================================================================================
FUNCTION read_nifti, fname, header=hdr, just_header=just_header,show=show, compress=compress, $
                            flip_ver=flip_ver, flip_hor=flip_hor


  ;; ------------------------------------------------------------------
  ;; check input
  ;; ------------------------------------------------------------------
  if (file_test(fname) ne 1) then return,-1
  if (strlen(fname)-3 eq strpos(fname,'.gz')) then compress=1  ; check for compressed file


  ;; ------------------------------------------------------------------
  ;; Create IDLffNifti Object
  ;; ------------------------------------------------------------------
  oNifti = obj_new('IDLffNifti')



  ;; ------------------------------------------------------------------
  ;; read header file
  ;; ------------------------------------------------------------------
  hdr = oNifti->ReadHeaderFile(fname, compress=compress)
  ;; show header ?
  if keyword_set(show) then begin
      print, '>>> filename:'
      print, fname
      oNifti->OutputHeader, hdr,-1
  endif
  ;; return ?
  if keyword_set(just_header) then begin
      obj_destroy, oNifti
      return, hdr
  endif


  ;; ------------------------------------------------------------------
  ;; read data and header
  ;; obj->ReadFile, [f] [,_EXTRA=extra])
  ;; Parameters:
  ;; fname       filename  (hdr, nii)
  ;; ------------------------------------------------------------------
  oNifti->ReadFile, fname, compress=compress
  vol = oNifti->GetData()
  ;; Note for IDL we need to flip top/bottom ==> this is the default 
  if not keyword_set(flip_ver) then vol = reverse(vol,2,/overwrite)

  ;; check if we need a left right flip (needed for HyperView)
  if keyword_set(flip_hor) then begin
     case hdr.dim[0] of
        2:  vol = reverse(vol,1,/overwrite)
        3:  vol = reverse(vol,3,/overwrite)
        else: stop
     endcase
  endif

  ;; ------------------------------------------------------------------
  ;; free pointer
  ;; ------------------------------------------------------------------
  obj_destroy, oNifti

  
  return, vol

END

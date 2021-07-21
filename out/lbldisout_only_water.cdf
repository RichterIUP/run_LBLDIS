CDF  �   
      n_wnum         n_cloud_layers        n_instances       n_height   E         Author        iLBLDIS developed by D.D. Turner, Space Science and Engineering Center, UW-Madison, dturner@ssec.wisc.edu       Version       G$Id: lblrtm_disort.c,v 1.44 2014/01/13 21:23:24 merrelli Release_3_0 $     Number_of_streams         16     RT_solver_number      0      Calculation_zenith_angle      180.000000 degrees      Calculation_zenith_angle_comment      I0 degrees is nadir (upwelling) while 180 degrees is zenith (downwelling)       
Input_flag        1      Path_to_LBLRTM_gas_ODs        3./out/.lblrtm_node25.cluster_20210716_14115252CEST     SSP_database_number_0         m/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/ssp_database/ssp_db.mie_wat.gamma_sigma_0p100       Wavenumber_comment        &Using selected microwindows option #0      Phase_function_comment        Using real phase functions     Solar_zenith_angle        2No solar source input included in the calculation      Solar_azimuth_angle       2No solar source input included in the calculation      Sun_earth_distance        2No solar source input included in the calculation      Solar_source_datafile         2No solar source input included in the calculation            wnum                	long_name         wavenumber     units         cm-1            �   radiance                   	long_name         Computed radiance      units         mW / (m2 sr cm-1)           �   flux_up                    	long_name         Computed upwelling flux    units         mW / (m2 cm-1)          �   	flux_down                      	long_name         "Computed downwelling diffuse flux      units         mW / (m2 cm-1)          �   	flux_beam                      	long_name         Computed direct beam flux      units         mW / (m2 cm-1)          �   
cld_height                 	long_name         Cloud height       units         km          p   	cld_dbnum                  	long_name         Cloud database number      units         	unitless       comment       0database number of SSP used for the cloud layer         t   cld_tau                   	long_name         Cloud optical depth    units         	unitless            x   cld_tau_wnum               	long_name         1Reference wavenumber for the cloud optical depth       units         cm-1       comment       ZA value of -1 implies the optical depth is the geometric limit value (i.e., where Qe = 2)           |   cld_reff               	long_name          Cloud particle effective radius    units         microns         �   sfc_temperature              	long_name         Surface temperature    units         K      comment       oA value of -1 indicates that the surface temperature is assumed to be the same as the lowest level temperature          �   sfc_emissivity                  	long_name         Surface emissivity     units         	unitless            �   pwv              	long_name         Precipitable water vapor       units         cm          �   height                 	long_name         Level height       units         km         �   pressure               	long_name         Level pressure     units         mb         �   temperature                	long_name         Level temperature      units         K          �   mixing_ratio               	long_name         Level water vapor mixing ratio     units         g/kg       comment       �This profile is estimated from the number of molecules per layer that is in the LBLRTM TAPE7 output.  This profile should be made to agree with the precipitable water vapor (PWV) amount that was also derived from the TAPE7, as the PWV is more accurate.           �=���  �    ��  @�  ��  ?�@<�t�=u=���>��>L��>�  >���>�33>�Q�>�p�>\>Ǯ>���>�ff?   ?�?�?�R?8Q�?G�?O\)?\(�?k�?���?�  ?�33?���?�ff@   @��@��@&ff@333@@  @L��@fff@s33@�  @�ff@���@�ff@�  @���@�33@���@�ff@�  @ٙ�@�33@���@�ffA   AffA��A33A��A   A(  A0  A8  A@  AH  AP  A`  Ap  A�  A�  A�  A�  D}hRD|{Dz��Dy@�Dw�HDv$Dt��Ds�Dr�Drr�Dr$�DqևDq��Dp�Dn��Dm�DlyyDk+Dh,JDfyHDe��Dd:�Db�D]�RDY��DU�HDPsdDKb=DFi�DA�yD<D8�D3��D/
-D*�bD"&VDD�PD9D.�D��D�D�9C�o�C�RC��CܤCӤZC��}CC�`bC���C�fFC��9C���C���C�G+Ct��Ca�CP�)CALC2�1C%��C[�B��FB�wLBi=qA��AQ��C���C���C�Y�C�&fC��3C�ٚC��3C�� C�s3C�s3C�s3C�s3C�s3C�&fC�� C�&fC�33C��3C�� C��3C��3C�� C��3C��fC���C�s3C��C�ffC��3C�  C�s3C��C���C�  C�Y�C�  C�Y�C���C�Y�C��fC��C�33C�33C~33C{�fCy��CwL�Cu�Cr�fCp��Cn�3ClffCi33Ce��Cb�C^��C[��CXL�CVL�CX��C\��C_� C`L�C`�fCaL�Cb� Cc33Cc�fCk�f@oK+@mN�@k��@iZ�@d΀@_�7@Z�@UK@Tj�@S�0@S�@Q�\@Qg@QQ�@Q��@Pk;@I0�@@�a@/��@$�V@{<@�Z@]?�Ul?�CY@��@(�@E��@J�#@9�{@$�@��@�@�z?���?��*?��2?c6-?&�U>��{>��v>�e�>���>�~�>��>��>��>��>�
>Y��>0	�>��=��U=���=��=�Z==G��=	�<�\�<~LE<��;��\;T�C;>�;5Z;/��;5|�;K�;X�
C�@ B���C��C��R    ?z�HC�� B�!�C�MC��    ?z�HC�@ B���C��	C��Z    ?z�HC�� B�G�C��.C���    ?z�HC�@ B���C�h�C�`    ?z�HC�� B�?"C���C��z    ?z�HC�@ B��$C�FhC�"    ?z�HC�� B��C��gC�{�    ?z�HC�@ B���C�"AC��2    ?z�HC�� B�I{C��.C�f    ?z�HC�@ B���C��`C��    ?z�HC�� B�=YC�iC�.^    ?z�HC�@ B���C��@C��~    ?z�HC�� B�V9C�AfC��    ?z�HC�@ B��C��WC���    ?z�HC�� B��bC�.C�>    ?z�HC�@ B��dC���C�Hu    ?z�HC�� B�>�C��C��J    ?z�HC�@ B��~C�TAC��t    ?z�HC�� B�I�C��SC�f     ?z�HC�@ B�sC�&zC���    ?z�HC�� B���C���C�e�    ?z�HC�@ B�P�C���C���    ?z�HC�� B��6C�^C�X    ?z�HC�@ B�l�C�ąC���    ?z�HC�� B��C�*sC��    ?z�HC�@ B�a�C��C�~    ?z�HC�� B���C���C��    ?z�HC�@ B�n;C�ZC�N    ?z�HC�� B��mC���C���    ?z�HC�@ B��C�!	C��    ?z�HC�� B�voC��C�:w    ?z�HC�@ B���C��C���    ?z�HC�� B��AC�I{C�    ?z�HC�@ B�V�C���C���    ?z�HC�� B��BC�fC� �    ?z�HC�@ B� C�m�C�@�    ?z�HC�� B�LKC��C�y    ?z�HC�@ B���C�,�C��    ?z�HC�� B�^XC���C�I�    ?z�HC�@ B�7#C��C�ۚ    ?z�HC�� B�jBC�J2C�    ?z�HC�@ B��ZC���C�l�    ?z�HC�� B�x3C��C��
    ?z�HC�@ B�b�C�a"C���    ?z�HC�� B�C��C�i�    ?z�HC�@ B�p�C��C���    ?z�HC�� B���C�uC�%�    ?z�HC�@ B��xC��JC���    ?z�HC�� B�b�C�,0C��    ?z�HC�@ B��C��C�w�    ?z�HC�� B�F�C�߀C���    ?z�HC�@ B���C�8NC�&j    ?z�HC�� B��C���C�}    ?z�HC�@ B£C���C��    ?z�HC�� B���C�?�C��    ?z�HC�@ B�P�C���C�v�    ?z�HC�� B���C��^C���    ?z�HC�@ B�@#C�CC�-�    ?z�HC�� B�-�C��$C�?@    ?z�HC�@ B�)[C��C�O(    ?z�HC�� B�}�C�>�C���    ?z�HC�@ B�/C��C� �    ?z�HC�� B�62C��C�    ?z�HC�@ Bŧ^C�7�C�Th    ?z�HC�� B��C���C��    ?z�HC�@ B�v�C�۸C���    ?z�HC�� B��C�-wC�]o    ?z�HC�@ B�{C�wC��5    ?z�HC�� BȰ�C�� C���    ?z�HC�@ B�� C�!�C��2    ?z�HC�� B��VC�o9C���    ?z�HC�@ B�+�C���C�:�    ?z�HC�� B���C�rC���    ?z�HC�@ B��C�\SC�N�    ?z�HC�� B�*�C��C��x    ?z�HC�@ B�p~C���C�ך    ?z�HC�� B�-6C�@TC�ʀ    ?z�HC�@ B˄C���C��    ?z�HC�� B�lC���C��c    ?z�HC�@ B�٘C��C�B�    ?z�HC�� B�.�C�i�C�z�    ?z�HC�@ ḂC���C��)    ?z�HC�� B͠�C��nC���    ?z�HC�@ B�b�C�IC�(�    ?z�HC�� B��qC��<C�v(    ?z�HC�@ B�B�C��FC�?V    ?z�HC�� B�P�C�rC�(�    ?z�HC�@ B���C�a�C�iV    ?z�HC�� B�C�C��0C��~    ?z�HC�@ BϪ�C��kC��w    ?z�HC�� B��C�2ZC�8�    ?z�HC�@ B� �C�w:C��B    ?z�HC�� B�C��#C�h�    ?z�HC�@ B���C��C���    ?z�HC�� B�ީC�D�C���    ?z�HC�@ BҚHC��8C�s    ?z�HC�� B���C��JC��{    ?z�HC�@ B�,+C�EC��*    ?z�HC�� B�Y�C�K�C��    ?z�HC�@ Bӛ�C���C�R�    ?z�HC�� B��C��>C��R    ?z�HC�@ Bӝ�C�C��    ?z�HC�� B�C!C�FkC�a    ?z�HC�@ BԬ�C��C���    ?z�HC�� BԚaC�¼C���    ?z�HC�@ B�9�C� �C�A�    ?z�HC�� B�a�C�=SC���    ?z�HC�@ BՁ�C�y|C���    ?z�HC�� B�C��DC��M    ?z�HC�@ B�.�C��lC�0�    ?z�HC�� B֣�C�+C�i    ?z�HC�@ B��C�eC���    ?z�HC�� B�PqC��HC��e    ?z�HC�@ B�yC���C�ޓ    ?z�HC�� B�
�C�yC��    ?z�HC�@ B׆�C�H�C��E    ?z�HC�� BأlC��.C�F�    ?z�HC�@ B�FFC��=C���    ?z�HC�� BَuC���C���    ?z�HC�@ Bل;C�'#C��0    ?z�HC�� B�GC�Z�C��5    ?z�HC�@ Bٔ*C���C��    ?z�HC�� B��$C���C�M�    ?z�HC�@ B�i5C��cC��V    ?z�HC�� B�qC�,�C��?    ?z�HC�@ B�\C�^�C���    ?z�HC�� B��C��C���    ?z�HC�@ B�&�C��C���    ?z�HC�� B۞�C��C��    ?z�HC�@ B�L�C�$�C�C�    ?z�HC�� B���C�V�C��1    ?z�HC�@ B�.�C��C���    ?z�HC�� B�@C���C��u    ?z�HC�@ B�}IC��]C�,u    ?z�HC�� B�0C�C� R    ?z�HC�@ B�5C�BxC��<    ?z�HC�� B�:C�l�C�q2    ?z�HC�@ B�l^C��	C���    ?z�HC�� B��YC��&C���    ?z�HC�@ BݬEC���C��`    ?z�HC�� B�(jC� �C��_    ?z�HC�@ B�j�C�K�C��    ?z�HC�� Bބ�C�vC��    ?z�HC�@ B�N�C��C��    ?z�HC�� B�R�C��yC��    ?z�HC�@ B��C��(C�i"    ?z�HC�� B�&jC��C�;A    ?z�HC�@ Bކ�C�AC�ku    ?z�HC�� B�%aC�h�C���    ?z�HC�@ B���C���C���    ?z�HC�� B��)C��C���    ?z�HC�@ B�@C��?C�;V    ?z�HC�� B��C��C�/�    ?z�HC�@ B��5C�%hC�k    ?z�HC�� B�tC�I�C�#5    ?z�HC�@ B��C�q6C��j    ?z�HC�� B��"C���C�w     ?z�HC�@ B�nC���C�`�    ?z�HC�� B��fC�ׅC��6    ?z�HC�@ B�O�C��C��Y    ?z�HC�� B�C�PC�p^    ?z�HC�@ B�#�C�=�C�b&    ?z�HC�� B��C�^C�z3    ?z�HC�@ B�#�C�z�C���    ?z�HC�� B��WC��3C�1X    ?z�HC�@ B��9C���C�ś    ?z�HC�� B��C��fC�`    ?z�HC�@ B��C���C�)    ?z�HC�� B�AC��C�/�    ?z�HC�@ B��RC�4oC�&    ?z�HC�� B��C�R�C��3    ?z�HC�@ B�RC�r�C�]y    ?z�HC�� B�nC���C��'    ?z�HC�@ B���C���C���    ?z�HC�� B�J�C���C��    ?z�HC�@ B�KC�ۮC��    ?z�HC�� B�XLC��oC��    ?z�HC�@ Bۯ�C��C�[b    ?z�HC�� B�`$C��C�V�    ?z�HC�@ B��C�>�C�     ?z�HC�� B�tC�XHC�e    ?z�HC�@ B��dC�p:C��:    ?z�HC�� B��C���C���    ?z�HC�@ B��`C��|C���    ?z�HC�� B�-�C��dC��S    ?z�HC�@ B�8C�ʵC���    ?z�HC�� B��C���C�Bt    ?z�HC�@ B�0C��(C��    ?z�HC�� Bݳ�C� #C�i    ?z�HC�@ B�0LC��C� �    ?z�HC�� B���C�1�C��    ?z�HC�@ B���C�DgC���    ?z�HC�� B��C�T�C� �    ?z�HC�@ B��|C�]C�N;    ?z�HC�� B�z�C�QC���    ?z�HC�@ B�c�C�^C��    ?z�HC�� B��eC�o�C��    ?z�HC�@ B��*C���C���    ?z�HC�� Bݣ�C���C��    ?z�HC�@ B�\C���C��;    ?z�HC�� B�H�C��!C���    ?z�HC�@ B��RC��C���    ?z�HC�� B�>C��eC��    ?z�HC�@ B��WC��C�u    ?z�HC�� B�'�C��C�4H    ?z�HC�@ B泞C� lC�j$    ?z�HC�� B�C�,=C�Z0    ?z�HC�@ B�qC�6�C�&5    ?z�HC�� B��C�<EC��e    ?z�HC�@ B��sC�N�C�t�    ?z�HC�� B�aWC�W)C�x    ?z�HD   B��IC�PC���    ?z�HD ` B�"C�UhC��    ?z�HD � B�&C�vC�~�    ?z�HD � B�_C�~�C��R    ?z�HD  B�H�C��C��g    ?z�HD` B�xC��5C� R    ?z�HD� B���C���C�a�    ?z�HD� B�xQC���C�I    ?z�HD  B�d�C���C���    ?z�HD` B�& C��oC��g    ?z�HD� B�M�C���C��:    ?z�HD� B��2C��C�i{    ?z�HD  B�!�C��:C���    ?z�HD` B�pC���C��    ?z�HD� B�vC��'C� ^    ?z�HD� B��C�ƞC�Y�    ?z�HD  B��C��C��Y    ?z�HD` B�}C���C��g    ?z�HD� B���C�s9C�O     ?z�HD� B�l�C�d�C�
L    ?z�HD  B��~C�^C��u    ?z�HD` B��NC���C���    ?z�HD� B��@C���C���    ?z�HD� B۾)C��WC��X    ?z�HD  B���C��BC��    ?z�HD` B�cC��!C�x�    ?z�HD� B�&C�r�C���    ?z�HD� B��C��?C�    ?z�HD  B���C��MC�Ա    ?z�HD` Bً�C���C���    ?z�HD� B�߅C�yC�i�    ?z�HD� B��C���C�ſ    ?z�HD  B�uC��jC���    ?z�HD` B�7C�ɻC���    ?z�HD� B现C�ƎC��+    ?z�HD� B�zC�C���    ?z�HD	  B��[C��QC�k    ?z�HD	` B⥻C��C��N    ?z�HD	� B�дC��C�[L    ?z�HD	� B��C�T�C�i�    ?z�HD
  B�-jC�f�C�$    ?z�HD
` B��~C�E�C���    ?z�HD
� B��C���C�/    ?z�HD
� B�P�C�ZC���    ?z�HD  B��C�
�C��}    ?z�HD` Bˈ�C�U�C���    ?z�HD� B�LC���C��    ?z�HD� B��tC���C�X�    ?z�HD  B��C�ѲC�Ҙ    ?z�HD` B��C��tC���    ?z�HD� B��C��C�A    ?z�HD� B�Z�C�C�C���    ?z�HD  B�iC��C��    ?z�HD` B�M�C���C���    ?z�HD� B��MC�$bC��}    ?z�HD� B��LC� �C�Y�    ?z�HD  B��C��C�h    ?z�HD` B��C��C�B�    ?z�HD� B�;�C��bC���    ?z�HD� Bۙ}C��XC� 6    ?z�HD  B�{PC�f�C�'    ?z�HD` B��C�O�C�D    ?z�HD� B�t�C���C���    ?z�HD� B�?C��:C��	    ?z�HD  B�v�C���C��o    ?z�HD` B�&�C�q�C��    ?z�HD� B���C�)6C��#    ?z�HD� B�^C�L�C�ܞ    ?z�HD  B���C�j1C��P    ?z�HD` B�P�C�R�C�:�    ?z�HD� B���C���C�R(    ?z�HD� B��C���C�7�    ?z�HD  B��]C�GC���    ?z�HD` B���C���C��R    ?z�HD� B��9C�e�C��-    ?z�HD� B�UC��C��>    ?z�HD  B�gC��C�/    ?z�HD` B��C�CC��T    ?z�HD� B��GC���C�p,    ?z�HD� B��C��jC��a    ?z�HD  B�/�C��C���    ?z�HD` B�xTC�OUC�z    ?z�HD� B�8C�T�C�|�    ?z�HD� Bθ�C�,MC��W    ?z�HD  B�C�C��C��    ?z�HD` B�&#C�I4C���    ?z�HD� B�l�C�H\C��t    ?z�HD� B�}�C��fC�x    ?z�HD  B�59C��C��&    ?z�HD` B�?vC���C���    ?z�HD� Bk�FC��aCtɵ    ?z�HD� B:[C���C�J    ?z�HD  B�TvC�b�C�o�    ?z�HD` Bk׼C�JpCt��    ?z�HD� BH��C��xCY��    ?z�HD� BN�(C��C^R�    ?z�HD  By�gC�wC~c�    ?z�HD` BM��C���C]Y�    ?z�HD� Ba��C���Cl��    ?z�HD� BT�kC��-Cbֺ    ?z�HD  B�scC�VGC��9    ?z�HD` B�9�C�9C�6    ?z�HD� B�BC��C��    ?z�HD� Bך7C�YC��8    ?z�HD  B�	�C�G�C���    ?z�HD` B�.�C���C�o    ?z�HD� Bp��C�/Cw�R    ?z�HD� B��QC�1�C�n�    ?z�HD  BЉ`C��JC�    ?z�HD` B��rC��C���    ?z�HD� BUa�C�\�Cb�S    ?z�HD� Bw�WC���C|o�    ?z�HD  B�}C�.�C��    ?z�HD` B�N�C�FC�3�    ?z�HD� B�@C���C�2z    ?z�HD� B]��C���Ci"�    ?z�HD  B:%�C�raCL�j    ?z�HD` B9�C�S.CL�Q    ?z�HD� B2#�C�#�CF#�    ?z�HD� B\j�C�[\Chq    ?z�HD  B��C�7C�VT    ?z�HD` B�	�C��C� �    ?z�HD� B�L�C���C�#e    ?z�HD� B�&C�
xC��    ?z�HD  B�#�C���C���    ?z�HD` B˝3C��4C�'S    ?z�HD� B�͹C�i�C�Տ    ?z�HD� B���C�3�C���    ?z�HD   B���C��~C��    ?z�HD ` BD��C��CU2�    ?z�HD � B�3�C�cVC���    ?z�HD � B���C�(�C��y    ?z�HD!  B��)C���C�cH    ?z�HD!` B�n�C�D	C�ђ    ?z�HD!� B(gC��bC=I�    ?z�HD!� B��C���C$pP    ?z�HD"  B
��C�p�C"g6    ?z�HD"` Bi�C�Q�C$�    ?z�HD"� B��C�>�C*��    ?z�HD"� B��C�a�C�CY    ?z�HD#  B$M7C�UC9��    ?z�HD#` B�]C���C��    ?z�HD#� A���C�z�C�    ?z�HD#� B��C�\�C�    ?z�HD$  B�C�DbC�    ?z�HD$` B�?C�W�C5<�    ?z�HD$� B�۱C�)	C���    ?z�HD$� B���C���C�)�    ?z�HD%  Br�2C��0Cv�    ?z�HD%` B(k�C���C<ˮ    ?z�HD%� B�J�C��EC�Z�    ?z�HD%� B���C�5�C���    ?z�HD&  BI�C�sC��    ?z�HD&` B )4C��-C    ?z�HD&� B%�C���C�-    ?z�HD&� BC�OC� �CSq    ?z�HD'  B�8�C���C��    ?z�HD'` B��5C���C�t�    ?z�HD'� Bt��C�LCwad    ?z�HD'� B���C�5�C� u    ?z�HD(  B!s4C�qC6[>    ?z�HD(` B�B~C��=C��    ?z�HD(� A��C�\XCK�    ?z�HD(� A�z�C� C��    ?z�HD)  A݂C��qC͓    ?z�HD)` B �;C��$C2�    ?z�HD)� B���C��C�{�    ?z�HD)� BB~C�/�CQz�    ?z�HD*  A�(�C�B�C�\    ?z�HD*` A�|�C��B�*    ?z�HD*� A���C��C[�    ?z�HD*� A��)C���B�Yl    ?z�HD+  A�6�C��TB���    ?z�HD+` A���C�r�C@    ?z�HD+� B�YC��CC1     ?z�HD+� B��GC���C��    ?z�HD,  Bģ�C��C�r1    ?z�HD,` B�~C�{�C�=�    ?z�HD,� B���C��C���    ?z�HD,� B��C���C"O�    ?z�HD-  Bk3C��yC%�h    ?z�HD-` B���C��]C��    ?z�HD-� B��aC�}~C�7N    ?z�HD-� BTz*C���C^��    ?z�HD.  B�4�C�I�C�Զ    ?z�HD.` B�F/C�"�C��    ?z�HD.� B��OC�vC�VB    ?z�HD.� A�@�C��Cx8    ?z�HD/  A�~C��4B��V    ?z�HD/` A�$�C���B���    ?z�HD/� A�8�C���Cp�    ?z�HD/� B��2C���C�_7    ?z�HD0  Bt^�C�,^Ct�b    ?z�HD0` B��C���C��$    ?z�HD0� Bb��C���ChT�    ?z�HD0� B��C�HC�q�    ?z�HD1  B�C��PC!��    ?z�HD1` A�*�C���B�F    ?z�HD1� A�8�C���B�d�    ?z�HD1� A��RC���B�'%    ?z�HD2  Aѭ�C��HC     ?z�HD2` B��C�nC���    ?z�HD2� B�C�|�C"}J    ?z�HD2� A��
C���B�H    ?z�HD3  A�xiC�q�B�9z    ?z�HD3` A�.�C�>�B�ï    ?z�HD3� B�BC��VC�-    ?z�HD3� A�!KC���B���    ?z�HD4  A��2C���B�A�    ?z�HD4` A�SC���C��    ?z�HD4� A�v�C�5�B�\n    ?z�HD4� A��C� 5B��~    ?z�HD5  A��C��TB�S�    ?z�HD5` A��eC��(B��u    ?z�HD5� A�b|C�n�B�L9    ?z�HD5� A���C�C�B�	�    ?z�HD6  A�UMC��ECT    ?z�HD6` B���C���C��O    ?z�HD6� B���C�AGCy�-    ?z�HD6� B4K�C��<CBԺ    ?z�HD7  A��bC�MDB���    ?z�HD7` A�/�C�� B���    ?z�HD7� A�(�C�˛B��    ?z�HD7� A{�3C���B�    ?z�HD8  Av6C�\�B��r    ?z�HD8` At)|C�'�B���    ?z�HD8� A�8rC���B��>    ?z�HD8� A�#�C���BҪ8    ?z�HD9  B�C�`]C��    ?z�HD9` A��jC���C ��    ?z�HD9� B_3C���Cb:�    ?z�HD9� A�>�C�l�Cz�    ?z�HD:  B�E�C��BC�qS    ?z�HD:` B�C�F�C�T    ?z�HD:� A�C�d B��    ?z�HD:� A�I�C�[�B�x\    ?z�HD;  Bs��C�|�Co�    ?z�HD;` A���C��WB�(�    ?z�HD;� Aj%�C�tlB���    ?z�HD;� A��iC�n]B���    ?z�HD<  B&��C�C6�O    ?z�HD<` B �"C���C1��    ?z�HD<� B�0MC�j�C���    ?z�HD<� A݃5C���C��    ?z�HD=  A�CC��B�B�    ?z�HD=` AN�LC��WB�@�    ?z�HD=� AH��C���B��    ?z�HD=� AE��C�s�B��p    ?z�HD>  A�d�C�v�B�k�    ?z�HD>` AU�:C�GB�[�    ?z�HD>� A> C�ʣB�A*    ?z�HD>� A;{�C���B}OH    ?z�HD?  A<1*C�[�B~!�    ?z�HD?` A<V�C�$�B~F    ?z�HD?� A<p�C���B~[�    ?z�HD?� A�.C��B���    ?z�HD@  A�+C�2�C��    ?z�HD@` A�X�C��C�    ?z�HD@� A�eC��C%H    ?z�HD@� A98�C���BzJq    ?z�HDA  A3�EC��JBs�7    ?z�HDA` A8+C�f�Bxۍ    ?z�HDA� AV�C�FB�K<    ?z�HDA� B>z'C�WVCG�    ?z�HDB  BQmC��cC,[4    ?z�HDB` BNC��rC"��    ?z�HDB� Ad��C�p�B���    ?z�HDB� BКC��C%>�    ?z�HDC  A9�C��Bz�:    ?z�HDC` A2V}C���Bq��    ?z�HDC� AB�C�u�B�9I    ?z�HDC� A��8C��B��    ?z�HDD  B���C��Cz�    ?z�HDD` A�#�C��B�ñ    ?z�HDD� A:B�C��0Bz��    ?z�HDD� A4#:C�O�Bs�I    ?z�HDE  A%r#C�HBa�t    ?z�HDE` A&-�C���Bbvi    ?z�HDE� A% �C��Ba#/    ?z�HDE� A=�mC�r�B~�k    ?z�HDF  A]<C�QB���    ?z�HDF` B�C��C�    ?z�HDF� A��*C�f�B�1M    ?z�HDF� BNC� yCPwc    ?z�HDG  B ��C�p�C.�d    ?z�HDG` B4GC��C!J�    ?z�HDG� B��/C��[C�"    ?z�HDG� BzNC��C��    ?z�HDH  A@�C�o�B���    ?z�HDH` A+	C�%�BhD    ?z�HDH� A��3C�Q�Bð�    ?z�HDH� BO��C�5�CQ|�    ?z�HDI  A3RC�~cBq��    ?z�HDI` A6��C�G4BvQ    ?z�HDI� AK�)C��B��X    ?z�HDI� AW��C��B�^Y    ?z�HDJ  B"zC���C$�J    ?z�HDJ` A�PC�B BI��    ?z�HDJ� A�gC��B;�@    ?z�HDJ� A�C��B9o>    ?z�HDK  A	^�C���B=��    ?z�HDK` A��#C�MB��	    ?z�HDK� B �!C�?C�    ?z�HDK� A�3C��BPO    ?z�HDL  AڣC��B=9�    ?z�HDL` @�?wC�^B0?�    ?z�HDL� @�(�C�"B.BG    ?z�HDL� @�	�C��B,�    ?z�HDM  @��C���B+�I    ?z�HDM` @�C�okB*�W    ?z�HDM� @�>4C�7�B/}v    ?z�HDM� @��QC���B*�j    ?z�HDN  A'�C���BH�a    ?z�HDN` A���C�B�[q    ?z�HDN� A
��C�U�B?jC    ?z�HDN� A��C��Cf    ?z�HDO  A�C��dBF�$    ?z�HDO` @�pC���B%<    ?z�HDO� @���C�X�B%��    ?z�HDO� @�	C��B"Zn    ?z�HDP  @��fC��\B"QI    ?z�HDP` @��C���B!��    ?z�HDP� @�]C�i�B"��    ?z�HDP� A(2C�W�BcEH    ?z�HDQ  @꒰C���B#�!    ?z�HDQ` @�{C��B#0X    ?z�HDQ� @��%C���B+r�    ?z�HDQ� A���C��B��    ?z�HDR  A�)/C�}yB�3�    ?z�HDR` A;q�C��By�-    ?z�HDR� A��C���BA    ?z�HDR� @���C�Y�B.�    ?z�HDS  @�fC��B0u    ?z�HDS` A �C���B6>�    ?z�HDS� A��C��HB3>�    ?z�HDS� A řC�m�B2�    ?z�HDT  Ae�C�8@B;�    ?z�HDT` A���C��C/_    ?z�HDT� A4C��}BLr8    ?z�HDT� A��C��/BX     ?z�HDU  BQ"�C��CM    ?z�HDU` A�PC�r<B��Y    ?z�HDU� A�˽C��dC�    ?z�HDU� A.�C���B8�d    ?z�HDV  @�rC�A�B�
    ?z�HDV` @��	C��B��    ?z�HDV� A�C���B>��    ?z�HDV� @���C��/B��    ?z�HDW  @���C�H�B�m    ?z�HDW` @�6�C��B��    ?z�HDW� @�*C��B��    ?z�HDW� @��*C&�B�    ?z�HDX  A?�C�BT��    ?z�HDX` A��iCnoB�N    ?z�HDX� @��OC}�qBĜ    ?z�HDX� @��C}C<B^�    ?z�HDY  @���C|��B,X    ?z�HDY` @��C|UB2�    ?z�HDY� @熗C|sB!�    ?z�HDY� AǲNC}dB�?�    ?z�HDZ  @�ACz�AB	��    ?z�HDZ` @��CzpB %�    ?z�HDZ� @���Cy�'A���    ?z�HDZ� @��Cy�NB�    ?z�HD[  @�y�Cy
B��    ?z�HD[` @���Cx�4B�    ?z�HD[� A�ոCyS�B��
    ?z�HD[� @��Cw�.B	c    ?z�HD\  @���Cw3B	c�    ?z�HD\` A	�1Cv��B<+<    ?z�HD\� @�F�Cv:�B
�    ?z�HD\� A�hCv�@B���    ?z�HD]  @ݣQCuhByg    ?z�HD]` @�(�CtˉA�U�    ?z�HD]� @�SCt��B$,v    ?z�HD]� A�S�Cu��B榦    ?z�HD^  A@N�Ct�B})�    ?z�HD^` @�pCr�*B�    ?z�HD^� @� /Cr�JB$��    ?z�HD^� ASCrDB6��    ?z�HD_  @��Cqv.A��    ?z�HD_` @�$MCp��A�F    ?z�HD_� @���Cp��A�t�    ?z�HD_� @���Cp�A�w    ?z�HD`  @�R>Co�A���    ?z�HD`` @�X�Co=�B��    ?z�HD`� @���Cn��A�"w    ?z�HD`� @�{BCn%�A��    ?z�HDa  @���Cm�mA��H    ?z�HDa` @�EACm4A��    ?z�HDa� @��$Cl�A�^�    ?z�HDa� @�f�ClC�Aԭ�    ?z�HDb  @�	:CǩA��    ?z�HDb` @��Ckh�A�ؔ    ?z�HDb� A%��CkoiB\�+    ?z�HDb� @�1�Cjp�A��	    ?z�HDc  A|�rCj��B�9�    ?z�HDc` AW��CjR<B�<J    ?z�HDc� @�GpCi1B��    ?z�HDc� @��%Ch�jA�    ?z�HDd  @��7Ch�A��    ?z�HDd` @��Cg��B{    ?z�HDd� @¾wCgGbB��    ?z�HDd� @���Cf�A�L    ?z�HDe  @�x�Cf7�Aݨ�    ?z�HDe` @�JCe�YA޵�    ?z�HDe� @�] Ce�HB(y     ?z�HDe� @�b+Cd�sA��    ?z�HDf  @���Cd]�A�ȑ    ?z�HDf` A<�CdQ,BC�    ?z�HDf� A�*Cd�DB���    ?z�HDf� @���Cc �A�΄    ?z�HDg  AWCb��BK!�    ?z�HDg` A6WnCb�7Boe    ?z�HDg� @��Ca��A���    ?z�HDg� @���Ca�AֺI    ?z�HDh  @���C`��A���    ?z�HDh` @���C`5UA�;�    ?z�HDh� @��C_��A��
    ?z�HDh� @�pSC_6eA�:j    ?z�HDi  @��%C^��Aɖ�    ?z�HDi` @�
C^>�A��|    ?z�HDi� @{�3C]�TA���    ?z�HDi� @t��C]J�A��R    ?z�HDj  @k�C\��A��    ?z�HDj` @r�:C\]�A�r    ?z�HDj� @a��C[�^A���    ?z�HDj� @^��C[iA��/    ?z�HDk  @�T�C[�A���    ?z�HDk` @��CZ��A�<�    ?z�HDk� @YOBCZ-A�	�    ?z�HDk� @\�dCY�A�h�    ?z�HDl  @�ѫCYF�A���    ?z�HDl` @a:�CX��A�d�    ?z�HDl� @�G�CXA�A�F�    ?z�HDl� @���CW�hA��    ?z�HDm  A�c�CX��B���    ?z�HDm` @cACV�6A��m    ?z�HDm� @N�#CVUA��    ?z�HDm� @K/CU�`A�s    ?z�HDn  @P��CUk�A�ں    ?z�HDn` @�&�CUo�B*F�    ?z�HDn� @wS�CT��A�.�    ?z�HDn� A��CT��BEit    ?z�HDo  @LECS��A���    ?z�HDo` @C�3CS;A���    ?z�HDo� @M�cCR�%A��e    ?z�HDo� @�>�CRg�A�]`    ?z�HDp  @U'CQȳA���    ?z�HDp` @U�^CQT�A�,m    ?z�HDp� @>{�CP�RA�=�    ?z�HDp� @:��CP`�A��=    ?z�HDq  @:L�CO�kA�R�    ?z�HDq` @r1CO��A�v�    ?z�HDq� @�PCO;A�3    ?z�HDq� @G��CN�	A��[    ?z�HDr  @6%iCN�A�aj    ?z�HDr` @5��CM��A��c    ?z�HDr� @d\CMH�A�Ώ    ?z�HDr� @��CM+
BI�    ?z�HDs  @G��CLU�A��E    ?z�HDs` @��CLMB�8    ?z�HDs� @���CK�JA�+�    ?z�HDs� @�?�CKl�B [�    ?z�HDt  @�CJ��B �    ?z�HDt` @�#vCJhB�o    ?z�HDt� @:*�CI��A� �    ?z�HDt� @*��CI%AzP�    ?z�HDu  @(�&CH��Aw��    ?z�HDu` @I"�CHM�A�=�    ?z�HDu� @9\�CGԠA�b    ?z�HDu� @,|�CG\�A|�    ?z�HDv  @[g}CF�bA��g    ?z�HDv` @#JCFt�Aoǌ    ?z�HDv� @!�CF�Al�-    ?z�HDv� @ 'wCE��Ak�7    ?z�HDw  @!�ECE�Amr(    ?z�HDw` @�CD�TAj
n    ?z�HDw� @-�CDA�A~��    ?z�HDw� @�CC�_Ae��    ?z�HDx  @��CCV�Ad}    ?z�HDx` @F�CB�cAcG�    ?z�HDx� @#��CBx�Ap�^    ?z�HDx� @{�CB�AfO    ?z�HDy  @�	CA��A_��    ?z�HDy` @+7CA#�Ag5�    ?z�HDy� @�#iC@�A�J-    ?z�HDy� @)H[C@H=Ax�    ?z�HDz  @�K�C@[A�j    ?z�HDz` @"D C?d�An8�    ?z�HDz� @BNC>�AY8�    ?z�HDz� @v��C>��A��)    ?z�HD{  @<j�C> �A�_    ?z�HD{` @�6C=�AV�    ?z�HD{� @C=-�AT�#    ?z�HD{� @,�C<��A{��    ?z�HD|  @��C<NAQO    ?z�HD|` @W��C;��A���    ?z�HD|� A1#C<P�BdS    ?z�HD|� @'�C;�Aux'    ?z�HD}  @
~,C:��AL�h    ?z�HD}` @��C:*$Af    ?z�HD}� A_>�C:�wB�vv    ?z�HD}� @��C9LlAc��    ?z�HD~  @K�C8�{AQ�    ?z�HD~` A�C9!�BC�4    ?z�HD~� @LHRC8�A���    ?z�HD~� @��C7�Aj�    ?z�HD  @?�C7�ABN&    ?z�HD` @\BC6��AW��    ?z�HD� @CVC6B:AG�    ?z�HD� ?�y&C5�4A<�W    ?z�HD� ?���C5c�A;�    ?z�HD�0 ?�H<C4�MA:��    ?z�HD�P ?�(�C4�}A;�    ?z�HD�p @]C4$+AS,:    ?z�HD�� Ai�C4f�B=�A    ?z�HD�� A{ C3�B2Q�    ?z�HD�� @�h0C3bA��E    ?z�HD�� @�rC2m�AE��    ?z�HD� @_�`C2(�A�a�    ?z�HD�0 ?���C1�{A6:    ?z�HD�P ?���C1#3A2�    ?z�HD�p ?�*�C0��A2�    ?z�HD�� ?�2�C0J�A/f    ?z�HD�� ?�A�C/��A74    ?z�HD�� @��C/��AW��    ?z�HD�� @�.C/ctA�t    ?z�HD� ?�|2C.�_A8"    ?z�HD�0 ?�"�C.2%A-�^    ?z�HD�P @q�TC-��A�?C    ?z�HD�p ?��C-[�A+}�    ?z�HD�� ?���C,�A)V�    ?z�HD�� ?��C,�A%��    ?z�HD�� ?��*C,�A%�    ?z�HD�� ?�\C+�A'h    ?z�HD� ?�vC+M�A<;    ?z�HD�0 @��C+P�B�    ?z�HD�P @��C*��AfYa    ?z�HD�p @�B�C*w!A�5�    ?z�HD�� ?�jC)�iA&��    ?z�HD�� ?�uC)6A#�Y    ?z�HD�� ?��;C(��A&K@    ?z�HD�� @��!C(�A�n�    ?z�HD� ?�C'��A(P�    ?z�HD�0 ?Ք�C'��A��    ?z�HD�P @-B	C'E�A{�o    ?z�HD�p @e�C&�AZ�~    ?z�HD�� @ʮ�C&��B	��    ?z�HD�� ?��C%�SA4@v    ?z�HD�� @ϟZC&qBs�    ?z�HD�� ?��C%'fA76    ?z�HD� ?�@�C$��A7[�    ?z�HD�0 @�pC$ǿB�    ?z�HD�P A�}�C%�RB��v    ?z�HD�p @(YC#��AS@    ?z�HD�� ?�C#lA'J�    ?z�HD�� ?�?gC"�A"�    ?z�HD�� ?���C"L�A�7    ?z�HD�� ?���C!��A m    ?z�HD� @�!aC!�&A�cI    ?z�HD�0 @L�C!/�Ac�    ?z�HD�P Am�4C!�PB��    ?z�HD�p @�C \AM�    ?z�HD�� @*xaC �Aw�    ?z�HD�� ?��C��Aa�    ?z�HD�� ?�(wC�A�    ?z�HD�� ?�,aC��A��    ?z�HD� ?�uCT*A'��    ?z�HD�0 @#C.Al�
    ?z�HD�P ?֛BC�`A>�    ?z�HD�p ?�:C$�A&�    ?z�HD�� ?�2C�A8��    ?z�HD�� @�a�C�dA�%�    ?z�HD�� ?��C�JA,�y    ?z�HD�� @�0C��AA9A    ?z�HD� @��C8�AD �    ?z�HD�0 ?�J�C�&A/�    ?z�HD�P @UFC�PA��J    ?z�HD�p Aj#C3?B�R�    ?z�HD�� @�C��A@i	    ?z�HD�� ?�N�C?�A.��    ?z�HD�� ?���CۆA+j    ?z�HD�� ?��rCx�A+�    ?z�HD� ?�Y{CRA3�    ?z�HD�0 ?���C��A6(Y    ?z�HD�P @�"C\�APM�    ?z�HD�p A_G�C�B��!    ?z�HD�� @aw^C��A��
    ?z�HD�� A�7�Cu�B��+    ?z�HD�� @;�C�A�z�    ?z�HD�� @�!Cl�A=�     ?z�HD� @dPCKAA�    ?z�HD�0 @!C�8Ae4    ?z�HD�P A���C�bB��    ?z�HD�p @�PC4�A�@    ?z�HD�� @CC��AWP�    ?z�HD�� ?�U�C#�A7�    ?z�HD�� @N��C�}A��    ?z�HD�� A��C�pB���    ?z�HD� @:^BCA���    ?z�HD�0 ?�fUC��A#p�    ?z�HD�P ?�m�C=�A)4    ?z�HD�p ?ɩCױA�    ?z�HD�� ?��ZCyA~Y    ?z�HD�� @Q`tCG�A��q    ?z�HD�� ?�۫C��A#��    ?z�HD�� ?��Ca�A-�f    ?z�HD� A7�C�WB`��    ?z�HD�0 A���C:.B�Sb    ?z�HD�P @	M�CM�AH:�    ?z�HD�p ?��C�A7BF    ?z�HD�� ?��C{�A�b    ?z�HD�� ?��C�A    ?z�HD�� ?�i<C��AC�    ?z�HD�� ?�5XC`)A f�    ?z�HD� ?��\CfA $�    ?z�HD�0 ?��qC��A!M�    ?z�HD�P ?�W�CKHA�    ?z�HD�p ?�SC
�A�^    ?z�HD�� ?�'�C
��A-_-    ?z�HD�� @%��C
S�AnA�    ?z�HD�� @���C
/�A���    ?z�HD�� A�A�C�.B���    ?z�HD� AS �C
2�Bz��    ?z�HD�0 A��C
e�B�@    ?z�HD�P @7w�C�sA���    ?z�HD�p ?�ISC�AZ    ?z�HD�� ?܈C��A"M�    ?z�HD�� @$�uCr^Am"�    ?z�HD�� ?�!4C��AE     ?z�HD�� ?�)�C�KA
�z    ?z�HD� ?��C=�@���    ?z�HD�0 ?��xC��A	)�    ?z�HD�P ?���C�SA�    ?z�HD�p ?�}C?hA7��    ?z�HD�� @0R�C��A{��    ?z�HD�� A��(C�B��    ?z�HD�� @Kc�CQ�A��<    ?z�HD�� AH�C�BnP!    ?z�HD� A�5C�AB��    ?z�HD�0 @,��C7�Av��    ?z�HD�P @e�CիA[�4    ?z�HD�p ?Ӊ�Cg�A�r    ?z�HD�� ?�%�C�A�F    ?z�HD�� ?�g�C��A�T    ?z�HD�� ?� CaxA(�n    ?z�HD�� ?ã�CAqu    ?z�HD� @cC �YAAzn    ?z�HD�0 @��C c�AM7I    ?z�HD�P ?�}C �A9��    ?z�HD�p @(�B�~�Aqr�    ?z�HD�� Ar8B��YB(!�    ?z�HD�� A�BC ��Bھ`    ?z�HD�� @��OB�<�B�G    ?z�HD�� @
��B���AH��    ?z�HD� ?���B�݆A�    ?z�HD�0 @9eB�txA�>G    ?z�HD�P @.�B���Aw�J    ?z�HD�p @OQ�B�,XA��    ?z�HD�� @��lB�� A��q    ?z�HD�� A�xjB��B��    ?z�HD�� B��B�G�B���    ?z�HD�� AAI(B�
�Bd��    ?z�HD� @�"gB�A�A��/    ?z�HD�0 @!�B��RAffD    ?z�HD�P @�r�B��&A�6    ?z�HD�p ?ߤ�B�}:A#�    ?z�HD�� @��B��_A��    ?z�HD�� @��B�>�AE`�    ?z�HD�� ?���B�zA��    ?z�HD�� @[OB���A<�    ?z�HD� @y�8B�hA�%H    ?z�HD�0 @���B�[A��j    ?z�HD�P A�*�B�B���    ?z�HD�p A�`TB�B��b    ?z�HD�� @�.�B�5rA��    ?z�HD�� @'gB�hAnG/    ?z�HD�� @UU�B핼A���    ?z�HD�� @�ɖB�e�A�W    ?z�HD� @WEB��A[}s    ?z�HD�0 @�G�B��-A��    ?z�HD�P ?�؝B�iA�    ?z�HD�p @(~B�3-Aoxa    ?z�HD�� @{�_B�қA�*    ?z�HD�� AE�B�rBf�    ?z�HD�� A��B�B��S    ?z�HD�� @��9B��rA��    ?z�HD� @T�B���AM�    ?z�HD�0 @MB�:`A<<�    ?z�HD�P ?�|@B�y�An�    ?z�HD�p ?�q�B���AA�    ?z�HD�� ?���B�<0AѺ    ?z�HD�� @\]B��9A���    ?z�HD�� @�PB��Bh�    ?z�HD�� @%UFB�Aj��    ?z�HD� @��B�?JA�    ?z�HD�0 @��BᢥA�    ?z�HD�P AdB���B3��    ?z�HD�p A���B��B���    ?z�HD�� A�C�B�=B��p    ?z�HD�� A	"�B�� B)�J    ?z�HD�� @�"B��'A��C    ?z�HD�� @܏@B�aXB H    ?z�HD� @P��B��A��`    ?z�HD�0 @�ͼB���A��K    ?z�HD�P Aޗ�B�lB��y    ?z�HD�p @�<�B��>B��    ?z�HD�� @��2B�i6B��    ?z�HD�� A#u�B�;�BC'�    ?z�HD�� @�9�B��BE�    ?z�HD�� @�*�B�L�A�Y�    ?z�HD� A��mB�ϹB��    ?z�HD�0 A�gZB��YB��    ?z�HD�P AA}�B�x�B^ƭ    ?z�HD�p @�:�B֕�A��    ?z�HD�� @%<`B՟Ai�    ?z�HD�� @���BղQA��    ?z�HD�� @ KSB�LEA8�    ?z�HD�� ?��0Bӡ�A��    ?z�HD� @/��B�BAv�    ?z�HD�0 @��Bғ�AU5�    ?z�HD�P ?�a�B��~A�    ?z�HD�p @�4B�U-A>��    ?z�HD�� @�փBѴ�B3�    ?z�HD�� A*��BѭBHD~    ?z�HD�� @Z��B��A�w�    ?z�HD�� Ai@B�	PB~`    ?z�HD� @�h�B�\�BX@    ?z�HD�0 @��zB�6�A�3�    ?z�HD�P ANIB�nBgu�    ?z�HD�p A���Bϸ!B�;R    ?z�HD�� A�{<BώB� :    ?z�HD�� A��B̯�B)m&    ?z�HD�� @�B�z�A�g�    ?z�HD�� @Je�Bʈ�A��X    ?z�HD� Ae�B���B$z    ?z�HD�0 @�/bB��,A�    ?z�HD�P @w��B��A�     ?z�HD�p A@r�BɷIBZ'�    ?z�HD�� AJMB�9�BbZ�    ?z�HD�� A�+B�D�B;��    ?z�HD�� A;�B���BU��    ?z�HD�� @�ŁB�V�A��    ?z�HD� @��B��GB �J    ?z�HD�0 A?�'B�J�BX�    ?z�HD�P A�O	B��lB�<�    ?z�HD�p Aڬ�BƼyB�*    ?z�HD�� A� �B�]B�|G    ?z�HD�� A���B�@�B���    ?z�HD�� A=OcB�t1BU��    ?z�HD�� AQ��B�`Bf�Y    ?z�HD� A~��B���B�ɾ    ?z�HD�0 A�.#B�Z�B��F    ?z�HD�P A��B��3B���    ?z�HD�p A���B���B��    ?z�HD�� A�dB��HB�jT    ?z�HD�� A�\B�3`B�JQ    ?z�HD�� A�owB��ZB��>    ?z�HD�� A�F)B�&#B���    ?z�HD� A��B���B��?    ?z�HD�0 A��/B�{vB�Ӕ    ?z�HD�P AɈB�dBB)i    ?z�HD�p @��9B��B�=    ?z�HD�� @��NB�B2�    ?z�HD�� @���B��A�S    ?z�HD�� @�gB���A�`�    ?z�HD�� Aq%B��NB|#F    ?z�HD� A�C�B���B��R    ?z�HD�0 A�]B��/B���    ?z�HD�P @�OB���B5k    ?z�HD�p ADB�BX��    ?z�HD�� A�t�B��MB�՚    ?z�HD�� A\n�B�(Bk�O    ?z�HD�� A���B���B�
    ?z�HD�� A�$jB�`�B�ٝ    ?z�HD� A�&IB���B���    ?z�HD�0 A�_aB�kB��]    ?z�HD�P A�I,B�}2B�J�    ?z�HD�p A4:�B���BJ^t    ?z�HD�� A5��B�'BKv�    ?z�HD�� AIB��@B[EL    ?z�HD�� A�B���B5>�    ?z�HD�� A>(:B���BQ��    ?z�HD� A�>�B�pB��C    ?z�HD�0 A�H7B�"�B�ia    ?z�HD�P @��B�r�B�h    ?z�HD�p @��B�ߺB	�{    ?z�HD�� @���B�k�Bf]    ?z�HD�� @��B�� A�U    ?z�HD�� @�I�B�^�B��    ?z�HD�� @�VbB���A�    ?z�HD� A
o�B��pB#	�    ?z�HD�0 A��B�biB�U    ?z�HD�P A�a�B���B�C�    ?z�HD�p A���B��B�t�    ?z�HD�� A�<GB�PNB�    ?z�HD�� A���B�\�B���    ?z�HD�� A�B��QB��6    ?z�HD�� A��{B��VB�(    ?z�HD� A��oB�N�B���    ?z�HD�0 A�g�B���B��c    ?z�HD�P A��8B�L�B��    ?z�HD�p AԴ�B���B�n�    ?z�HD�� AՍ2B�J]B��Z    ?z�HD�� A�t�B��
B��j    ?z�HD�� AҚ�B�H�B��    ?z�HD�� A�U�B���B�O    ?z�HD� A��B�H�B�q�    ?z�HD�0 A��3B��B��(    ?z�HD�P A��B�B��q    ?z�HD�p A�;
B���B�7�    ?z�HD�� A��VB�fB�'�    ?z�HD�� A���B�)iB��J    ?z�HD�� A��B��FB�@�    ?z�HD�� A�6B�z`B�l    ?z�HD� A��~B��^B��    ?z�HD�0 A��B��B�i?    ?z�HD�P A���B�pB�˴    ?z�HD�p A|��B�1*By    ?z�HD�� A��_B�X�B���    ?z�HD�� A�V�B��TB���    ?z�HD�� A��GB�U�B��S    ?z�HD�� A�R&B��B��'    ?z�HD� A���B�|9B�(y    ?z�HD�0 A��zB���B�Y    ?z�HD�P A���B��B���    ?z�HD�p Aƺ�B��B�k    ?z�HD�� A�D�B���B��    ?z�HD�� A�4�B��B�H    ?z�HD�� A��B��bB�r�    ?z�HD�� A�bB�#QB���    ?z�HD� A���B���B���    ?z�HD�0 A�9B��B�o�    ?z�HD�P A��RB�yLB�+�    ?z�HD�p A��/B�=�B�`    ?z�HD�� A�ĽB���B�O    ?z�HD�� A���B�\!B���    ?z�HD�� A�D�B��B�'�    ?z�HD�� A�\B�W!B�$u    ?z�HD� A�SB��oB���    ?z�HD�0 A�>�B�t�B�`    ?z�HD�P A��+B�B���    ?z�HD�p A�eDB���B�    ?z�HD�� A�\nB�)CB�h�    ?z�HD�� A� HB���B���    ?z�HD�� A���B�<KB���    ?z�HD�� A��B�ʂB��7    ?z�HD� A�B�S�B��    ?z�HD�0 A��lB��B�͡    ?z�HD�P A� �B�s�B��u    ?z�HD�p A���B�7B�X#    ?z�HD�� A���B��sB�@A    ?z�HD�� A�B��B�9�    ?z�HD�� A��#B��fB�u�    ?z�HD�� A��B�4B�z�    ?z�HD� A�?�B���B���    ?z�HD�0 A�pB�OB�a    ?z�HD�P A��B��B�i�    ?z�HD�p A���B�n�B�J�    ?z�HD�� A���B� �B�?k    ?z�HD�� A��B��KB��    ?z�HD�� A���B�(.B��    ?z�HD�� A��LB��VB�߭    ?z�HD� A��XB�@4B��    ?z�HD�0 A��B���B���    ?z�HD�P A�/B�goB���    ?z�HD�p A���B��ZB�T{    ?z�HD�� A�,�B���B���    ?z�HD�� A�B�#�B��"    ?z�HD�� A��B��'B��u    ?z�HD�� A�ssB�D�B��p    ?z�HD� A���B��B���    ?z�HD�0 A���B�g�B��    ?z�HD�P A���B�=B�ZC    ?z�HD�p A�B��gB�N�    ?z�HD�� A���B�)�B��    ?z�HD�� A�j�B���B�Qh    ?z�HD�� A��;B�PVB�f�    ?z�HD�� A��B��oB��    ?z�HD� A�wkB�y$B�$y    ?z�HD�0 A�WLB�kB���    ?z�HD�P A��HB��fB�~�    ?z�HD�p A���B�FB��+    ?z�HD�� A��B�ٍB��*    ?z�HD�� A��B�qwB���    ?z�HD�� A�ƖB�
B�|-    ?z�HD�� A��B���B�N�    ?z�HD� A�6�B�4�B��    ?z�HD�0 A�>B��$B~�t    ?z�HD�P A�ބB�`XB}�N    ?z�HD�p A��B��@B~�    ?z�HD�� A�X�B���B��    ?z�HD�� A��6B�4B    ?z�HD�� A���B�RB{&�    ?z�HD�� A���B~��Bz�7    ?z�HD� A�L�B}�By�s    ?z�HD�0 A���B}-=By�H    ?z�HD�P A��B|`Bx�I    ?z�HD�p A�6�B{�5ByRI    ?z�HD�� A�igBz�\BwY�    ?z�HD�� A���Bz �Bu�    ?z�HD�� A��By6~Bt�    ?z�HD�� A�F�Bxs�Bu�    ?z�HD� A�rBw��Bv�    ?z�HD�0 A���Bv�Bv8�    ?z�HD�P A��Bv,"Buv�    ?z�HD�p A��Buh�Bu8    ?z�HD�� A���Bt��Bs��    ?z�HD�� A���Bs��BqD�    ?z�HD�� A�WBsxBp�    ?z�HD�� A�[}BrQ�Bqd�    ?z�HD� A�o�Bq�HBp/    ?z�HD�0 A�1Bp��Bl��    ?z�HD�P A���Bo�"Bk��    ?z�HD�p A��Bo6
Bj�:    ?z�HD�� A���Bn�VBlv�    ?z�HD�� A�Y�BmĽBl�|    ?z�HD�� A��BmBk��    ?z�HD�� A�9Bl<�Bi��    ?z�HD� A�l�Bk|�Bh�/    ?z�HD�0 A��Bj��Bid�    ?z�HD�P A���Bj[Bg�W    ?z�HD�p A���BiKABhV    ?z�HD�� A��'Bh��BhI'    ?z�HD�� A��Bg�nBgA;    ?z�HD�� A�^%Bg�Be:5    ?z�HD�� A��xBfMcBc �    ?z�HD� A��gBe�yBa�y    ?z�HD�0 A���BdӧB`é    ?z�HD�P A�E�BdwB_��    ?z�HD�p A��MBc^`B_�    ?z�HD�� A�WsBb��B^f�    ?z�HD�� A�	&Ba�]B^��    ?z�HD�� A��|BaA�B_s�    ?z�HD�� A�+�B`��B_��    ?z�HD� A��8B_�B^�Q    ?z�HD�0 A��B_GB\p    ?z�HD�P A���B^`�B[m    ?z�HD�p A�AoB]�1B[�|    ?z�HD�� A�F�B\��B[�    ?z�HD�� A�k�B\F�BZT�    ?z�HD�� A��B[�=BZ��    ?z�HD�� A�`BZ�BZ�    ?z�HD� A�� BZ5*BZ    ?z�HD�0 A�O�BY��BYX    ?z�HD�P A��%BX�JBX�-    ?z�HD�p A�NBX$BW��    ?z�HD�� A��:BWg(BU�W    ?z�HD�� A��#BV�8BT�    ?z�HD�� A�]
BVQBS�n    ?z�HD�� A��oBUW�BT    ?z�HD� A�BT�BTJ�    ?z�HD�0 A�x�BS��BS|�    ?z�HD�P A�C�BSKWBQ��    ?z�HD�p A�WJBR��BP��    ?z�HD�� A�rBQ��BO��    ?z�HD�� A�dBQ@QBOI\    ?z�HD�� A��BP��BO�I    ?z�HD�� A�BO�BO�    ?z�HD� A���BOF*BO�    ?z�HD�0 A�`kBN��BNr�    ?z�HD�P A�v�BM�BMEU    ?z�HD�p A���BMD`BL�_    ?z�HD�� A��BL��BL1=    ?z�HD�� A�UJBK��BJYT    ?z�HD�� A~`BK?�BH��    ?z�HD�� A|�`BJ�PBG��    ?z�HD� A}�BI�HBH�    ?z�HD�0 A~֢BIR�BH��    ?z�HD�P Az�BH�mBF%F    ?z�HD�p Ay��BG�iBE�    ?z�HD�� AyT�BGW�BD��    ?z�HD�� Ay��BF�rBD�    ?z�HD�� A{'�BFBE�>    ?z�HD�� Az�QBEu*BE*]    ?z�HD� AyU�BDφBD(<    ?z�HD�0 AyxfBD.�BD�    ?z�HD�P Ax��BC��BCY@    ?z�HD�p AvۄBB��BB:�    ?z�HD�� At��BB@�B@��    ?z�HD�� As�CBA�B@+    ?z�HD�� As��BA B?��    ?z�HD�� At�B@cB?��    ?z�HD� As��B?�wB?��    ?z�HD�0 ArkXB?#HB>�B    ?z�HD�P Aq��B>�B>7�    ?z�HD�p Ap�<B=�PB=b�    ?z�H
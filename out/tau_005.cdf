CDF     
      n_wnum         n_cloud_layers     
   n_instances       n_height   E         Author        iLBLDIS developed by D.D. Turner, Space Science and Engineering Center, UW-Madison, dturner@ssec.wisc.edu       Version       G$Id: lblrtm_disort.c,v 1.44 2014/01/13 21:23:24 merrelli Release_3_0 $     Number_of_streams         16     RT_solver_number      0      Calculation_zenith_angle      180.000000 degrees      Calculation_zenith_angle_comment      I0 degrees is nadir (upwelling) while 180 degrees is zenith (downwelling)       
Input_flag        1      Path_to_LBLRTM_gas_ODs        {/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/run_LBLDIS/out/.lblrtm_node52.cluster_20210813_16204848CEST     SSP_database_number_0         m/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/ssp_database/ssp_db.mie_wat.gamma_sigma_0p100       SSP_database_number_1         m/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/ssp_database/ssp_db.mie_ice.gamma_sigma_0p100       Wavenumber_comment        &Using selected microwindows option #0      Phase_function_comment        Using real phase functions     Solar_zenith_angle        2No solar source input included in the calculation      Solar_azimuth_angle       2No solar source input included in the calculation      Sun_earth_distance        2No solar source input included in the calculation      Solar_source_datafile         2No solar source input included in the calculation            wnum                	long_name         wavenumber     units         cm-1            X   radiance                   	long_name         Computed radiance      units         mW / (m2 sr cm-1)           \   flux_up                    	long_name         Computed upwelling flux    units         mW / (m2 cm-1)          `   	flux_down                      	long_name         "Computed downwelling diffuse flux      units         mW / (m2 cm-1)          d   	flux_beam                      	long_name         Computed direct beam flux      units         mW / (m2 cm-1)          h   
cld_height                 	long_name         Cloud height       units         km        (  L   	cld_dbnum                  	long_name         Cloud database number      units         	unitless       comment       0database number of SSP used for the cloud layer         t   cld_tau                   	long_name         Cloud optical depth    units         	unitless          (  �   cld_tau_wnum               	long_name         1Reference wavenumber for the cloud optical depth       units         cm-1       comment       ZA value of -1 implies the optical depth is the geometric limit value (i.e., where Qe = 2)         (  �   cld_reff               	long_name          Cloud particle effective radius    units         microns       (  �   sfc_temperature              	long_name         Surface temperature    units         K      comment       oA value of -1 indicates that the surface temperature is assumed to be the same as the lowest level temperature              sfc_emissivity                  	long_name         Surface emissivity     units         	unitless            l   pwv              	long_name         Precipitable water vapor       units         cm             height                 	long_name         Level height       units         km            pressure               	long_name         Level pressure     units         mb            temperature                	long_name         Level temperature      units         K          0   mixing_ratio               	long_name         Level water vapor mixing ratio     units         g/kg       comment       �This profile is estimated from the number of molecules per layer that is in the LBLRTM TAPE7 output.  This profile should be made to agree with the precipitable water vapor (PWV) amount that was also derived from the TAPE7, as the PWV is more accurate.           D=���=���>L��>L��>���>���>���>���?   ?                  =���    =���    =���    =���    =���    ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  @�  A�  @�  A�  @�  A�  @�  A�  @�  A�  ��      <�t�=u=���>��>L��>�  >���>�33>�Q�>�p�>\>Ǯ>���>�ff?   ?�?�?�R?8Q�?G�?O\)?\(�?k�?���?�  ?�33?���?�ff@   @��@��@&ff@333@@  @L��@fff@s33@�  @�ff@���@�ff@�  @���@�33@���@�ff@�  @ٙ�@�33@���@�ffA   AffA��A33A��A   A(  A0  A8  A@  AH  AP  A`  Ap  A�  A�  A�  A�  D}hRD|{Dz��Dy@�Dw�HDv$Dt��Ds�Dr�Drr�Dr$�DqևDq��Dp�Dn��Dm�DlyyDk+Dh,JDfyHDe��Dd:�Db�D]�RDY��DU�HDPsdDKb=DFi�DA�yD<D8�D3��D/
-D*�bD"&VDD�PD9D.�D��D�D�9C�o�C�RC��CܤCӤZC��}CC�`bC���C�fFC��9C���C���C�G+Ct��Ca�CP�)CALC2�1C%��C[�B��FB�wLBi=qA��AQ��C���C���C�Y�C�&fC��3C�ٚC��3C�� C�s3C�s3C�s3C�s3C�s3C�&fC�� C�&fC�33C��3C�� C��3C��3C�� C��3C��fC���C�s3C��C�ffC��3C�  C�s3C��C���C�  C�Y�C�  C�Y�C���C�Y�C��fC��C�33C�33C~33C{�fCy��CwL�Cu�Cr�fCp��Cn�3ClffCi33Ce��Cb�C^��C[��CXL�CVL�CX��C\��C_� C`L�C`�fCaL�Cb� Cc33Cc�fCk�f                                                                                                                                                                                                                                                                                    |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  C� BQEC��hCq;    ?z�HC� B�UC���C�v    ?z�HC� B�C�I:C     ?z�HC� B,�C���Ch\    ?z�HC� Bs�C��_C��    ?z�HC� B��C��C �    ?z�HC� B�C�[�C Qu    ?z�HC� B@�C���C ��    ?z�HC� B�C�قC �=    ?z�HC� B�xC��C!+�    ?z�HC� B �C�O�C!qk    ?z�HC� B=�C��C!��    ?z�HC� B{�C��oC!��    ?z�HC� B�C��C"9�    ?z�HC� B�C�%�C"x�    ?z�HC�� B)JC�V�C"��    ?z�HC� Ba�C���C"�    ?z�HC� B��C���C#.o    ?z�HC� B�KC��C#hT    ?z�HC� B �C��C#�5    ?z�HC�� B3}C�/?C#��    ?z�HC�� BeyC�U*C$
�    ?z�HC�� B��C�yJC$>N    ?z�HC�� B��C���C$o�    ?z�HC�� B�KC��*C$�8    ?z�HC�� B	�C���C$§    ?z�HC�� B	#BC���C$ֻ    ?z�HC�� B	2�C��C$�G    ?z�HC�� B	AsC�,C$��    ?z�HC�� B	OfC�C�C%
o    ?z�HC�� B	ZjC�YnC%0    ?z�HD @ B	f�C�m�C%%Z    ?z�HD � B	qC��C%0b    ?z�HD@ B	zPC���C%;�    ?z�HD� B	��C���C%C�    ?z�HD@ B	�4C��;C%J    ?z�HD� B	�!C���C%P�    ?z�HD@ B	��C��C%U    ?z�HD� B	��C��xC%W�    ?z�HD@ B	�{C��KC%Y�    ?z�HD� B	��C��uC%Y�    ?z�HD@ B	��C��C%X~    ?z�HD� B	�3C���C%V    ?z�HD@ B	�bC��gC%Q�    ?z�HD� B	��C��4C%Lj    ?z�HD@ B	�C��hC%E|    ?z�HD� B	��C��C%=�    ?z�HD@ B	~�C��(C%4)    ?z�HD� B	w~C�ǼC%)�    ?z�HD	@ B	oC���C%c    ?z�HD	� B	ZVC��2C%"    ?z�HD
@ B	:pC���C$�    ?z�HD
� B	�C��RC$�Z    ?z�HD@ B�C��C$�H    ?z�HD� B�3C�z}C$uE    ?z�HD@ B� C�h^C$M�    ?z�HD� B��C�T�C$$�    ?z�HD@ Be�C�?�C#�g    ?z�HD� B=�C�)XC#�c    ?z�HD@ B�C�}C#�    ?z�HD� B��C��.C#u�    ?z�HD@ B��C��C#Gj    ?z�HD� B�6C��nC#�    ?z�HD@ BlfC���C"�}    ?z�HD� B@4C��-C"��    ?z�HD@ B�C�e C"�H    ?z�HD� B�aC�C|C"L�    ?z�HD@ B��C� �C"�    ?z�HD� B��C��sC!�}    ?z�HD@ BR�C���C!�Q    ?z�HD� B �C��,C!pK    ?z�HD@ B��C��C!7T    ?z�HD� B�TC�^�C ��    ?z�HD@ B�4C�4+C ��    ?z�HD� BPC�FC ��    ?z�HD@ B�C��C =<    ?z�HD� BʸC���C�    ?z�HD@ B��C�|�C��    ?z�HD� B:&C�K�CP�    ?z�HD@ B�C��C�J    ?z�HD� B��C��{C�%    ?z�HD@ B\�C���C[    ?z�HD� B�C�|RC    ?z�HD@ B�lC�ExC�v    ?z�HD� BwRC��C^    ?z�HD@ B+'C��kC�    ?z�HD� BܾC��6C��    ?z�HD@ B��C�^�CX5    ?z�HD� B=oC�"vC�l    ?z�HD@ B �*C���C�:    ?z�HD� B ��C��\CK    ?z�HD@ B K{C�f�C�v    ?z�HD� A���C�%�C��    ?z�HD@ A�N7C��CC6�    ?z�HD� A���C��}C�W    ?z�HD @ A��nC�]�Cz	    ?z�HD � A�VYC��CE    ?z�HD!@ A���C��C��    ?z�HD!� A� )C��=CZ�    ?z�HD"@ A�S�C�DxC�t    ?z�HD"� A���C���C�g    ?z�HD#@ A���C���Cs    ?z�HD#� A���C�gC�    ?z�HD$@ A�.�C�MC7�    ?z�HD$� A�_5C�ήC�    ?z�HD%@ A���C��CO�    ?z�HD%� A���C�2�C۩    ?z�HD&@ A��C��JCe    ?z�HD&� A��C��C�#    ?z�HD'@ A�B�C�A�Cx6    ?z�HD'� A�mJC��
C �    ?z�HD(@ A�C��<C��    ?z�HD(� A�C�I�CP    ?z�HD)@ A��PC��-C�b    ?z�HD)� A��C���C�    ?z�HD*@ A�6�C�I�C�    ?z�HD*� A�ZEC���C'1    ?z�HD+@ A�MC��JC�[    ?z�HD+� A��C�B�C/�    ?z�HD,@ A�ƮC��C�)    ?z�HD,� A��^C���C6~    ?z�HD-@ A�
�C�4�C��    ?z�HD-� A�+�C�هC:�    ?z�HD.@ A�LfC�}kC�q    ?z�HD.� A�i�C� �C>@    ?z�HD/@ A�C���C��    ?z�HD/� A�mC�d�C1�    ?z�HD0@ A�+C��C
��    ?z�HD0� A�JC���C
<    ?z�HD1@ A�C�E�C	�    ?z�HD1� A༁C��C	�    ?z�HD2@ A��;C���C�d    ?z�HD2� A�̏C� �C��    ?z�HD3@ A�ڏC���Ck�    ?z�HD3� A��nC�Z7C�@    ?z�HD4@ A��AC��CTE    ?z�HD4� A��)C��[C��    ?z�HD5@ A��NC�,C;�    ?z�HD5� A� C��'C��    ?z�HD6@ A�"C�_�C"4    ?z�HD6� A�?C���C�;    ?z�HD7@ A��C��7C|    ?z�HD7� A�0C�)=Cy�    ?z�HD8@ A�!OC���C�    ?z�HD8� A�#�C�WsC]�    ?z�HD9@ A�*~C���C ��    ?z�HD9� A�/�C���C @Q    ?z�HD:@ A�5C�B�b�    ?z�HD:� A�:qC��B�D�    ?z�HD;@ A�=IC�BuB�&�    ?z�HD;� A�*0C��XB��    ?z�HD<@ A��C�i�B���    ?z�HD<� A��xC��vB�W�    ?z�HD=@ AɴmC���B�U    ?z�HD=� AȏC� �B��G    ?z�HD>@ A�eC��rB�o4    ?z�HD>� A�:�C�C�B� �    ?z�HD?@ A��C��cB��    ?z�HD?� A��C�d�B�    ?z�HD@@ A¿qC���B�3�    ?z�HD@� A���C��AB��    ?z�HDA@ A�k�C�zB�F    ?z�HDA� A�B�C��bB�A�    ?z�HDB@ A�<C�0�B��    ?z�HDB� A��[C���B鉼    ?z�HDC@ A���C�L�B�-�    ?z�HDC� A�mzC��(B��1    ?z�HDD@ A�:mC�gQB�u�    ?z�HDD� A�2C��]B�F�    ?z�HDE@ A�~�C���B�iv    ?z�HDE� A�z�C�=B�8@    ?z�HDF@ A���C���B�<�    ?z�HDF� A�$C�)qB�    ?z�HDG@ A�c�C��%B��    ?z�HDG� A��{C�=�B���    ?z�HDH@ A���C��Bۆ�    ?z�HDH� A���C�P�Bلl    ?z�HDI@ A�%RC���B��Z    ?z�HDI� A���C�d�B�3�    ?z�HDJ@ A�lVC��hB�    ?z�HDJ� A�&�C�xB�Q    ?z�HDK@ A���C�dB��N    ?z�HDK� A���C���B�q�    ?z�HDL@ A�YC��B��    ?z�HDL� A�)C��^B͔t    ?z�HDM@ A��]C�$�B�%�    ?z�HDM� A��GC��dBʸ    ?z�HDN@ A�M�C�5�B�J�    ?z�HDN� A���C���B�º    ?z�HDO@ A���C�E~B�A    ?z�HDO� A��\C�� B�a�    ?z�HDP@ A��C�U�B�6�    ?z�HDP� A��C��B�7�    ?z�HDQ@ A�g�C�f4B�co    ?z�HDQ� A���C��B�Ρ    ?z�HDR@ A�J�C�w�B���    ?z�HDR� A��UC�B��    ?z�HDS@ A���C���B��*    ?z�HDS� A���C��B�,�    ?z�HDT@ A�#�C���B�,!    ?z�HDT� A���C��B���    ?z�HDU@ A��FC��AB�<�    ?z�HDU� A���C�#�B�T    ?z�HDV@ A�'LC��B�V�    ?z�HDV� A��C�-B��/    ?z�HDW@ A���C���B�Q�    ?z�HDW� A�@C�7�B��S    ?z�HDX@ A��C{�B�
�    ?z�HDX� A�q�C~�B���    ?z�HDY@ A�[�C}��B�[�    ?z�HDY� A�~C|��B��;    ?z�HDZ@ A��C{�
B�uC    ?z�HDZ� A}@�Cz�DB�>    ?z�HD[@ A{�Cy��B�ѯ    ?z�HD[� Ax��Cx��B�\G    ?z�HD\@ Av�Cw�B���    ?z�HD\� At%�Cv�B��)    ?z�HD]@ Ap�	Cu��B��9    ?z�HD]� An5CuB���    ?z�HD^@ Al��Ct�B��    ?z�HD^� Ai5tCs�B��V    ?z�HD_@ Af"8Cr&)B��r    ?z�HD_� Ac(8Cq1�B�@Y    ?z�HD`@ A`G^Cp=�B��+    ?z�HD`� A]t�CoI�B�    ?z�HDa@ AZ�CnV�B���    ?z�HDa� AY#Cmd^B���    ?z�HDb@ AW��Clr�B��v    ?z�HDb� AV�Ck��B��    ?z�HDc@ AV^�Cj�.B��    ?z�HDc� AV��Ci�_B���    ?z�HDd@ AV��Ch�:B���    ?z�HDd� AUR?CgÚB�F    ?z�HDe@ ARדCf�}B��s    ?z�HDe� AQ�0Ce�dB�    ?z�HDf@ AUVCd��B��h    ?z�HDf� AR�`CdB�~.    ?z�HDg@ AI'eCc>B�7�    ?z�HDg� AJ8&CbB��S    ?z�HDh@ AI�Ca0�B�Z�    ?z�HDh� AF"�C`>�B�m�    ?z�HDi@ A@S�C_I�B�=c    ?z�HDi� A:��C^U;BzG:    ?z�HDj@ A5��C]a�Bt��    ?z�HDj� A3C\q�Bq��    ?z�HDk@ A1!WC[�XBo�7    ?z�HDk� A/A�CZ�Bmb�    ?z�HDl@ A-`�CY�Bk:_    ?z�HDl� A+|TCX�WBi�    ?z�HDm@ A)��CW��Bf�    ?z�HDm� A(#�CVߎBe9P    ?z�HDn@ A'%�CU�BdG    ?z�HDn� A&hCU�Bb��    ?z�HDo@ A%�CTBa��    ?z�HDo� A$CS3`B`p0    ?z�HDp@ A"�ZCRI#B_.    ?z�HDp� A!վCQ_8B]�    ?z�HDq@ A ��CPu�B\�    ?z�HDq� A�hCO�QB[D�    ?z�HDr@ Ah@CN�cBY�    ?z�HDr� A5�CM��BX��    ?z�HDs@ AZCLҙBW)O    ?z�HDs� A��CK�BU�    ?z�HDt@ A��CK@BTQ�    ?z�HDt� AN�CJ)BR��    ?z�HDu@ A�CI5cBQf9    ?z�HDu� A�+CHOBO�3    ?z�HDv@ Ag�CGiBNd    ?z�HDv� A�CF��BLڸ    ?z�HDw@ A�RCE�cBKNz    ?z�HDw� AZ�CD��BI��    ?z�HDx@ A��CC�ZBH&'    ?z�HDx� A�NCB�vBF��    ?z�HDy@ A)�CBBD�    ?z�HDy� A
�CA+BCC�    ?z�HDz@ A	�C@H�BB    ?z�HDz� A	�C?g�BA=f    ?z�HD{@ Ac"C>�fB@j�    ?z�HD{� A��C=�AB?��    ?z�HD|@ AC<ǮB>��    ?z�HD|� A[0C;�B=�V    ?z�HD}@ A�FC;	�B=r    ?z�HD}� A�C:+�B<$w    ?z�HD~@ A8HC9M�B;?`    ?z�HD~� A{�C8p�B:X�    ?z�HD@ A�C7��B9o�    ?z�HD� A]C6��B8�C    ?z�HD�  A>�C5��B7�8    ?z�HD�` A zC5 �B6��    ?z�HD�� @�a�C4&B5�b    ?z�HD�� @���C3K�B4�    ?z�HD�  @�=#C2rFB3�    ?z�HD�` @���C1�$B2�    ?z�HD�� @��MC0�xB1�7    ?z�HD�� @�Y�C/�zB0��    ?z�HD�  @��tC/�B/��    ?z�HD�` @��C.9�B.��    ?z�HD�� @�RvC-c�B-�    ?z�HD�� @�C,��B,��    ?z�HD�  @��C+�@B+�y    ?z�HD�` @�LC*��B*�f    ?z�HD�� @�:C*B*>    ?z�HD�� @�CC)<�B)e�    ?z�HD�  @�jC(jWB(�,    ?z�HD�` @��C'�WB(
�    ?z�HD�� @�4�C&�eB'�D    ?z�HD�� @�1C%��B'A)    ?z�HD�  @�C%&�B&��    ?z�HD�` @�"C$W�B&~�    ?z�HD�� @�XdC#�!B&-�    ?z�HD�� @�֔C"�B%�U    ?z�HD�  @�,�C!�dB%E$    ?z�HD�` @�PDC! =B$��    ?z�HD�� @��C S�B$�    ?z�HD�� @�TwC�@B#�F    ?z�HD�  @�C�dB$gM    ?z�HD�` @�;�C� B$�    ?z�HD�� @�4C,B%#0    ?z�HD�� @��FCa�B$�    ?z�HD�  @�/C�/B#�3    ?z�HD�` @��C�+B#�    ?z�HD�� @�SC	�B"֣    ?z�HD�� @�UWCC�B"�?    ?z�HD�  @��ZC~B"�/    ?z�HD�` @�\�C��B"f�    ?z�HD�� @�C�B!�    ?z�HD�� @�C0�B"ql    ?z�HD�  @��Ch{B��    ?z�HD�` @���C�uBKY    ?z�HD�� @��oC�Bxf    ?z�HD�� @ػgC�B:�    ?z�HD�  @խ�C\Bpq    ?z�HD�` @�b�C��B�    ?z�HD�� @ϫ~C�`B�    ?z�HD�� @͊CB�E    ?z�HD�  @�*�CYB�    ?z�HD�` @��C��B�T    ?z�HD�� @��@C�*B9�    ?z�HD�� @�A�C wB��    ?z�HD�  @ǭ�Cc�B�'    ?z�HD�` @��)C�UB�    ?z�HD�� @ŇC
�3B5�    ?z�HD�� @�C
2�B�>    ?z�HD�  @Ä�C	yJB�    ?z�HD�` @��C��BI(    ?z�HD�� @�c�C�B��    ?z�HD�� @�N�CQ�B�    ?z�HD�  @¡JC�~B�+    ?z�HD�` @�pC�B¬    ?z�HD�� @��C1)B��    ?z�HD�� @�|C|hB%�    ?z�HD�  @�h�C��Bi    ?z�HD�` @��CdB
ȸ    ?z�HD�� @��8CbB
��    ?z�HD�� @�C��B��    ?z�HD�  @�l,C�B��    ?z�HD�` @�.]C KTBIs    ?z�HD�� @�1�B�5�B��    ?z�HD�� @�b�B���B"    ?z�HD�  @�}�B�x�B^�    ?z�HD�` @�(�B��B��    ?z�HD�� @�{�B��FBt    ?z�HD�� @��$B�d*B]=    ?z�HD�  @�mzB�	�B ��    ?z�HD�` @���B��MA�؇    ?z�HD�� @���B�[lA���    ?z�HD�� @�IB��A�w%    ?z�HD�  @�B0B�A�G=    ?z�HD�` @�s�B�a�A��    ?z�HD�� @��RB�^A���    ?z�HD�� @�� B��%A���    ?z�HD�  @��B�tWA��%    ?z�HD�` @�"sB�'�A�T�    ?z�HD�� @�S�B�ܖA�&    ?z�HD�� @��XB蒻A��"    ?z�HD�  @��B�JEA��*    ?z�HD�` @�	�B�:A���    ?z�HD�� @�O�B�~A�     ?z�HD�� @��VB�yA�    ?z�HD�  @��B�6A�    ?z�HD�` @�5�B��RA�q�    ?z�HD�� @�{�B߳�A�]�    ?z�HD�� @��B�t�A�J5    ?z�HD�  @�B�7
A�6�    ?z�HD�` @�\�B���A�&�    ?z�HD�� @��>Bڿ�A�c    ?z�HD�� @��Bم�A��    ?z�HD�  @�>�B�MPA���    ?z�HD�` @���B�+A��    ?z�HD�� @��B��dA��`    ?z�HD�� @� �Bԫ�A��    ?z�HD�  @�n�B�x�A�Q    ?z�HD�` @��EB�F�A�    ?z�HD�� @��B�8Aߟu    ?z�HD�� @�WB��Aޓi    ?z�HD�  @���BιA݇%    ?z�HD�` @��B͌dA�z�    ?z�HD�� @�A�B�a2A�q     ?z�HD�� @��aB�7$A�dY    ?z�HD�  @��B��A�\�    ?z�HD�` @�1qB���A�T�    ?z�HD�� @���B���A�J�    ?z�HD�� @��XBƜ3A�Dp    ?z�HD�  @�+�B�x�A�@�    ?z�HD�` @��jB�V}A�9�    ?z�HD�� @��{B�5�A�5�    ?z�HD�� @�.B�A�1�    ?z�HD�  @���B���A�.�    ?z�HD�` @�ܙB�ڦA�-�    ?z�HD�� @�0�B���A�+m    ?z�HD�� @��.B��~A�-    ?z�HD�  @��RB��IA�*�    ?z�HD�` @�8�B�sqA�+�    ?z�HD�� @��~B�\�A�,�    ?z�HD�� @��B�G�A�.�    ?z�HD�  @�C�B�3eA�1u    ?z�HD�` @��yB� �A�4    ?z�HD�� @��B�A�:.    ?z�HD�� @�S�B���A�>�    ?z�HD�  @���B���A�E�    ?z�HD�` @�0B��A�I�    ?z�HD�� @�gEB�ՙA�R~    ?z�HD�� @���B��bA�Zx    ?z�HD�  @�jB��aA�c�    ?z�HD�` @�z�B���A�m    ?z�HD�� @���B��%A�{/    ?z�HD�� @�?�B���A��`    ?z�HD�  @��DB���A���    ?z�HD�` @�
�B��4A���    ?z�HD�� @�o{B���A��i    ?z�HD�� @���B��gA���    ?z�HD�  @�6�B��aA��J    ?z�HD�` @���B���A��
    ?z�HD�� @��B���A�k    ?z�HD�� @~��B���A�"    ?z�HD�  @}��B��mA�0�    ?z�HD�` @|s�B��kA�F�    ?z�HD�� @{F�B���A�]    ?z�HD�� @z�B��4A�wN    ?z�HD�  @x�B���A��k    ?z�HD�` @w��B���A��    ?z�HD�� @v�SB���A��t    ?z�HD�� @uh,B�� A��    ?z�HD�  @t>�B���A��    ?z�HD�` @s7B�/A��    ?z�HD�� @q�cB�!A�67    ?z�HD�� @p��B�."A�T�    ?z�HD�  @o�fB�AfA�s�    ?z�HD�` @n{�B�U�A��`    ?z�HD�� @mU�B�kHA���    ?z�HD�� @l1�B���A���    ?z�HD�  @k%B���A���    ?z�HD�` @i�^B���A�V    ?z�HD�� @h�)B��'A�@�    ?z�HD�� @g�B��A�e�    ?z�HD�  @f�TB�A���    ?z�HD�` @eiB�"�A���    ?z�HD�� @dH�B�A�A���    ?z�HD�� @c),B�a~A���    ?z�HD�  @b�B���A�*�    ?z�HD�` @`�B���A�T�    ?z�HD�� @_�B��?A�}<    ?z�HD�� @^ŚB��A��    ?z�HD�  @]�B�NA��o    ?z�HD�` @\�jB�9%A��    ?z�HD�� @[�B�`�A�//    ?z�HD�� @Zl�B���A�]�    ?z�HD�  @YS�B��A��    ?z�HD�` @XHYB��0A��s    ?z�HD�� @W/�B��A���    ?z�HD�� @V�B�8�A��    ?z�HD�  @UxB�gaA�Q�    ?z�HD�` @TB���A���    ?z�HD�� @R�BB�ǬA��    ?z�HD�� @Q�DB��fA��@    ?z�HD�  @P��B�,PA�,P    ?z�HD�` @O�{B~��A�r�    ?z�HD�� @OB}*�A��     ?z�HD�� @N/B{��A��    ?z�HD�  @MIBzxA�K    ?z�HD�` @Lb�Bxu�A��(    ?z�HD�� @KxBv�A�ܲ    ?z�HD�� @J�Bu]hA�&�    ?z�HD�  @I��Bs�7A�q�    ?z�HD�` @H�	BrMA��l    ?z�HD�� @G�Bp�7A�y    ?z�HD�� @F��BoE\A�U$    ?z�HD�  @F�BmāA���    ?z�HD�` @E9BlE�A��}    ?z�HD�� @DXBj�A�>5    ?z�HD�� @Cr$BiNwA��a    ?z�HD�  @B�BgջA���    ?z�HD�` @A��Bf_1A�)    ?z�HD�� @@�1Bd��A�z�    ?z�HD�� @?�QBcx&A��>    ?z�HD�  @?{Bb�A�h    ?z�HD�` @>6B`�A�p�    ?z�HD�� @=T�B_,�A��    ?z�HD�� @<w�B]�A��    ?z�HD�  @;��B\Y�A�g�    ?z�HD�` @:��BZ��A��2    ?z�HD�� @9�BY�cA��    ?z�HD�� @9�BX+�A�h�    ?z�HD�  @8+�BV��A���    ?z�HD�` @7LZBUlA�~    ?z�HD�� @6t�BTOA~��    ?z�HD�� @5��BR�mA}��    ?z�HD�  @4�BQ[{A|;�    ?z�HD�` @3�`BPPAz�    ?z�HD�� @3yBN�Ay��    ?z�HD�� @2=>BM[�AxX�    ?z�HD�  @1d�BL
�Aw�    ?z�HD�` @0��BJ� Auċ    ?z�HD�� @/�
BImWAt�    ?z�HD�� @.�EBH!�As9�    ?z�HD�  @.�BF׶Aq�;    ?z�HD�` @->3BE��Ap��    ?z�HD�� @,j�BDIiAom�    ?z�HD�� @+��BC�An/(    ?z�HD�  @*�$BA�uAl��    ?z�HD�` @)�YB@��Ak��    ?z�HD�� @)+)B?B�AjqS    ?z�HD�� @([
B>�Ai5,    ?z�HD�  @'��B<�NAg��    ?z�HD�` @&��B;��Af��    ?z�HD�� @&!9B:Y$Ae��    ?z�HD�� @%��B9#�Ad�<    ?z�HD�  @%J\B7��Ad�    ?z�HD�` @$�YB6��Ac1�    ?z�HD�� @$uB5��Ab\C    ?z�HD�� @$�B4_OAa�q    ?z�HD�  @#�VB32nA`�4    ?z�HD�` @#&�B2WA_Ѐ    ?z�HD�� @"�^B0�
A^�    ?z�HD�� @"B�B/�jA^�    ?z�HD�  @!�4B.�gA]A    ?z�HD�` @!^�B-l	A\e�    ?z�HD�� @ �PB,I�A[��    ?z�HD�� @ w�B+(�AZ�Z    ?z�HD�  @ �B*	uAY�O    ?z�HD�` @��B(��AX�R    ?z�HD�� @�B'��AX�    ?z�HD�� @�LB&��AW@�    ?z�HD�  @$�B%��AVe    ?z�HD�` @�cB$��AU�R    ?z�HD�� @1�B#p�AT��    ?z�HD�� @��B"\�AS�[    ?z�HD�  @>(B!J�AR�    ?z�HD�` @��B :#AR�    ?z�HD�� @E:B+!AQ8]    ?z�HD�� @R9BhAP�_    ?z�HD�  @�MB�AP��    ?z�HD�` @jUB
�AQQ    ?z�HD  @�6B@AQ�    ?z�HD�� @`B�[AQ$    ?z�HD�  @��B��AQ&�    ?z�HD�` @;"B�
AQ�    ?z�HDà @�sB�AQ6    ?z�HD�� @�B��AQ�    ?z�HD�  @T�B�8AP��    ?z�HD�` @��B�BAP�    ?z�HDĠ @��B��AP��    ?z�HD�� @ 2�B�AP|�    ?z�HD�  @ qB;APN<    ?z�HD�` @ ��BAP>    ?z�HDŠ @ ִBwAO��    ?z�HD�� @!�B*JAO�c    ?z�HD�  @!)oB7�AOX/    ?z�HD�` @!N�BFNAO[    ?z�HDƠ @!joBV\AN��    ?z�HD�� @!%B
g�ANo{    ?z�HD�  @!��B	z�AN�    ?z�HD�` @!��B�<AM��    ?z�HDǠ @!��B��AMZM    ?z�HD�� @!��B�AL�F    ?z�HD�  @!��BԏAL��    ?z�HD�` @!��B�XAK�h    ?z�HDȠ @!k�B	�AKrI    ?z�HD�� @!D%B%�AJ�r    ?z�HD�  @! BC�AJL�    ?z�HD�` @ ��BcAI��    ?z�HDɠ @ ��B ��AI�    ?z�HD�� @ �A�KAHqY    ?z�HD�  @ CA���AGʺ    ?z�HD�` @�EA���AG7    ?z�HDʠ @��A�&zAFqK    ?z�HD�� @hyA�t�AE�B    ?z�HD�  @�A���AE�    ?z�HD�` @�%A��ADKS    ?z�HDˠ @a>A�o�AC��    ?z�HD�� @L�A���AB*�    ?z�HD�  @��A�A@,
    ?z�HD�` @��A�y�A>-�    ?z�HD̠ @_A��7A<6$    ?z�HD�� @X{A�7)A:D�    ?z�HD�  @�3A��A8U�    ?z�HD�` @�A���A6lE    ?z�HD͠ @O	A�f�A4��    ?z�HD�� @�WA��A2��    ?z�HD�  @�A�=�A0ɒ    ?z�HD�` @lFA�/A.�    ?z�HDΠ @
ՏA�A-/    ?z�HD�� @	@�Aޓ�A+P    ?z�HD�  @��A�
XA)��    ?z�HD�` @%}Aۃ�A'��    ?z�HDϠ @��A��NA%�!    ?z�HD�� @A�}{A$<�    ?z�HD�  @��A��cA"��    ?z�HD�` @ �AՁfA ̕    ?z�HDР ?�N�A��A�    ?z�HD�� ?�T�AҎ�Am    ?z�HD�  ?�M�A��A��    ?z�HD�` ?�FZAϥ{A�    ?z�HDѠ ?�AfA�4|Arp    ?z�HD�� ?�;�A�ŧAɝ    ?z�HD�  ?�AqA�YA$�    ?z�HD�` ?�@�A���A}    ?z�HDҠ ?�?�AȆ�A�+    ?z�HD�� ?�HZA�!FA6�    ?z�HD�  ?�Q\AŽ�A��    ?z�HD�` ?�X�A�\�A��    ?z�HDӠ ?�]�A���AP�    ?z�HD�� ?�klA��A	��    ?z�HD�  ?ӄA�F�A    ?z�HD�` ?Џ�A��Au    ?z�HDԠ ?ͧQA���A��    ?z�HD�� ?ʸ�A�C�A;�    ?z�HD�  ?��SA��:A�?    ?z�HD�` ?��A���A \    ?z�HDՠ ?�	�A�T�@�޺    ?z�HD�� ?�,
A�	h@��-    ?z�HD�  ?�L:A��@���    ?z�HD�` ?�r)A�x�@�[�    ?z�HD֠ ?���A�3�@�2'    ?z�HD�� ?���A��@�/    ?z�HD�  ?���A���@��    ?z�HD�` ?�NA�p�@��0    ?z�HDנ ?�C�A�3�@��    ?z�HD�� ?�x(A���@���    ?z�HD�  ?��+A���@�n�    ?z�HD�` ?��A��@�V�    ?z�HDؠ ?�-A�T@�=�    ?z�HD�� ?�U-A�!/@�&i    ?z�HD�  ?���A��Z@�    ?z�HD�` ?���A��K@�#    ?z�HD٠ ?��A��b@��G    ?z�HD�� ?�f,A�iT@��?    ?z�HD�  ?���A�@;@��    ?z�HD�` ?���A��@���    ?z�HDڠ ?�B�A��@���    ?z�HD�� ?��<A��!@��    ?z�HD�  ?��#A���@��7    ?z�HD�` ?�9fA���@�ʫ    ?z�HD۠ ?%AA�q3@�Β    ?z�HD�� ?y��A�U+@��    ?z�HD�  ?t�QA�;@�֦    ?z�HD�` ?oKPA�"�@���    ?z�HDܠ ?j!A�s@��    ?z�HD�� ?d��A���@���    ?z�HD�  ?_�4A���@��Z    ?z�HD�` ?ZtjA���@�h    ?z�HDݠ ?UN�A�ā@��    ?z�HD�� ?P%?A��@�.    ?z�HD�  ?KA��A@�A^    ?z�HD�` ?E�4A��P@�W\    ?z�HDޠ ?@�0A��@�nI    ?z�HD�� ?;��A��|@���    ?z�HD�  ?6�bA���@���    ?z�HD�` ?1�OA���@���    ?z�HDߠ ?,�!A��9@��t    ?z�HD�� ?'��A��S@��    ?z�HD�  ?"��A��X@|)    ?z�HD�` ?�#A���@vi�    ?z�HD� ?�]A��$@p�9    ?z�HD�� ?ϋA��@j��    ?z�HD�  ?ҰA���@g<[    ?z�HD�` ?�bA���@e��    ?z�HD� ?��A��@c��    ?z�HD�� ?�QA���@b(�    ?z�HD�  ?�A��@`~�    ?z�HD�` ?sA�`@^�<    ?z�HD� ?
p6A~�@]/�    ?z�HD�� ?	c�A|?�@[��    ?z�HD�  ?ZAzoX@Y��    ?z�HD�` ?V�Ax�Q@XS     ?z�HD� ?T�Av�)@V��    ?z�HD�� ?M,Au)@U!2    ?z�HD�  ?L�AsMP@S��    ?z�HD�` ?K�Aq�h@Q�    ?z�HD� ?N�Ao�t@Ph�    ?z�HD�� ?S�An�@Nߍ    ?z�HD�  ? ]1Al[�@MV�    ?z�HD�` >��2Aj��@Kʽ    ?z�HD� >��Ah��@JH�    ?z�HD�� >��QAgE�@H�_    ?z�HD�  >���Ae�u@GE    ?z�HD�` >�*?Ac�-@E��    ?z�HD� >�>tAbI�@DN�    ?z�HD�� >�jA`�L@Bݠ    ?z�HD�  >�A_�@Ad    ?z�HD�` >�VA]g�@?�_    ?z�HD� >��tA[�{@>��    ?z�HD�� >��AZ4L@=�    ?z�HD�  >�D�AX��@;�v    ?z�HD�` >�~�AW�@:B    ?z�HD� >淊AU{�@8��    ?z�HD�� >���AS�@7x    ?z�HD�  >�-^ARc�@6�    ?z�HD�` >�AP�@4��    ?z�HD� >�ǇAOV�@3_�    ?z�HD�� >��AM� @2s    ?z�HD�  >�TjALT@0�4    ?z�HD�` >ڪfAJ��@/YS    ?z�HD� >��9AI[�@.�    ?z�HD�� >�Q�AG��@,�(    ?z�HD�  >ղAFn)@+j�    ?z�HD�` >��AD�O@*�    ?z�HD� >�d�AC��@(��    ?z�HD�� >�ŵAB�@'��    ?z�HD�  >�*qA@�B@&G+    ?z�HD�` >͔1A?HX@%�    ?z�HD� >���A=��@#Ř    ?z�HD�� >�lWA<}�@"�?    ?z�HD�  >��A;@!L    ?z�HD�` >�O
A9�9@ {    ?z�HD�� >�0�A8`�@U    ?z�HD�� >�x�A7�@K�    ?z�HD�  >��VA5��@�-    ?z�HD�` >�A4ZV@��    ?z�HD� >�YA3�@��    ?z�HD�� >tA1�T@/�    ?z�HD�  >��IA0i,@h�    ?z�HD�` >�&�A/�@�1    ?z�HD� >�t�A-�@�    ?z�HD�� >���A,��@g    ?z�HD�  >��,A+G�@U[    ?z�HD�` >�C�A*Y@��    ?z�HD� >��A(��@ϧ    ?z�HD�� >���A'��@    ?z�HD�  >�A&J�@Ca    ?z�HD�` >�OjA%�@�    ?z�HD� >��A#�4@�{    ?z�HD�� >���A"��@��    ?z�HD�  >�%A!pm@9    ?z�HD�` >�SA ?H@w�    ?z�HD� >���A%@�    ?z�HD�� >��^A�.@�K    ?z�HD�  >��A�[@7�    ?z�HD�` >�P,A��@wG    ?z�HD� >���Ah�@�C    ?z�H
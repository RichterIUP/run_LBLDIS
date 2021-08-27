CDF     
      n_wnum         n_cloud_layers     
   n_instances       n_height   E         Author        iLBLDIS developed by D.D. Turner, Space Science and Engineering Center, UW-Madison, dturner@ssec.wisc.edu       Version       G$Id: lblrtm_disort.c,v 1.44 2014/01/13 21:23:24 merrelli Release_3_0 $     Number_of_streams         16     RT_solver_number      0      Calculation_zenith_angle      180.000000 degrees      Calculation_zenith_angle_comment      I0 degrees is nadir (upwelling) while 180 degrees is zenith (downwelling)       
Input_flag        1      Path_to_LBLRTM_gas_ODs        {/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/run_LBLDIS/out/.lblrtm_node52.cluster_20210813_16155454CEST     SSP_database_number_0         m/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/ssp_database/ssp_db.mie_wat.gamma_sigma_0p100       SSP_database_number_1         m/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/ssp_database/ssp_db.mie_ice.gamma_sigma_0p100       Wavenumber_comment        &Using selected microwindows option #0      Phase_function_comment        Using real phase functions     Solar_zenith_angle        2No solar source input included in the calculation      Solar_azimuth_angle       2No solar source input included in the calculation      Sun_earth_distance        2No solar source input included in the calculation      Solar_source_datafile         2No solar source input included in the calculation            wnum                	long_name         wavenumber     units         cm-1            X   radiance                   	long_name         Computed radiance      units         mW / (m2 sr cm-1)           \   flux_up                    	long_name         Computed upwelling flux    units         mW / (m2 cm-1)          `   	flux_down                      	long_name         "Computed downwelling diffuse flux      units         mW / (m2 cm-1)          d   	flux_beam                      	long_name         Computed direct beam flux      units         mW / (m2 cm-1)          h   
cld_height                 	long_name         Cloud height       units         km        (  L   	cld_dbnum                  	long_name         Cloud database number      units         	unitless       comment       0database number of SSP used for the cloud layer         t   cld_tau                   	long_name         Cloud optical depth    units         	unitless          (  �   cld_tau_wnum               	long_name         1Reference wavenumber for the cloud optical depth       units         cm-1       comment       ZA value of -1 implies the optical depth is the geometric limit value (i.e., where Qe = 2)         (  �   cld_reff               	long_name          Cloud particle effective radius    units         microns       (  �   sfc_temperature              	long_name         Surface temperature    units         K      comment       oA value of -1 indicates that the surface temperature is assumed to be the same as the lowest level temperature              sfc_emissivity                  	long_name         Surface emissivity     units         	unitless            l   pwv              	long_name         Precipitable water vapor       units         cm             height                 	long_name         Level height       units         km            pressure               	long_name         Level pressure     units         mb            temperature                	long_name         Level temperature      units         K          0   mixing_ratio               	long_name         Level water vapor mixing ratio     units         g/kg       comment       �This profile is estimated from the number of molecules per layer that is in the LBLRTM TAPE7 output.  This profile should be made to agree with the precipitable water vapor (PWV) amount that was also derived from the TAPE7, as the PWV is more accurate.           D=���=���>L��>L��>���>���>���>���?   ?                  ?���    ?���    ?���    ?���    ?���    ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  @�  A�  @�  A�  @�  A�  @�  A�  @�  A�  ��      <�t�=u=���>��>L��>�  >���>�33>�Q�>�p�>\>Ǯ>���>�ff?   ?�?�?�R?8Q�?G�?O\)?\(�?k�?���?�  ?�33?���?�ff@   @��@��@&ff@333@@  @L��@fff@s33@�  @�ff@���@�ff@�  @���@�33@���@�ff@�  @ٙ�@�33@���@�ffA   AffA��A33A��A   A(  A0  A8  A@  AH  AP  A`  Ap  A�  A�  A�  A�  D}hRD|{Dz��Dy@�Dw�HDv$Dt��Ds�Dr�Drr�Dr$�DqևDq��Dp�Dn��Dm�DlyyDk+Dh,JDfyHDe��Dd:�Db�D]�RDY��DU�HDPsdDKb=DFi�DA�yD<D8�D3��D/
-D*�bD"&VDD�PD9D.�D��D�D�9C�o�C�RC��CܤCӤZC��}CC�`bC���C�fFC��9C���C���C�G+Ct��Ca�CP�)CALC2�1C%��C[�B��FB�wLBi=qA��AQ��C���C���C�Y�C�&fC��3C�ٚC��3C�� C�s3C�s3C�s3C�s3C�s3C�&fC�� C�&fC�33C��3C�� C��3C��3C�� C��3C��fC���C�s3C��C�ffC��3C�  C�s3C��C���C�  C�Y�C�  C�Y�C���C�Y�C��fC��C�33C�33C~33C{�fCy��CwL�Cu�Cr�fCp��Cn�3ClffCi33Ce��Cb�C^��C[��CXL�CVL�CX��C\��C_� C`L�C`�fCaL�Cb� Cc33Cc�fCk�f                                                                                                                                                                                                                                                                                    |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  C� B�̓C��C��p    ?z�HC� B�+zC���C��    ?z�HC� B݆�C�3C�8#    ?z�HC� B��C�{;C�~�    ?z�HC� B�6�C��tC��    ?z�HC� Bފ�C�C��    ?z�HC� B���C�G�C�Fd    ?z�HC� B�,4C��GC��    ?z�HC� B�y=C�ƟC��$    ?z�HC� B���C�2C��    ?z�HC� B�@C�=�C�6I    ?z�HC� B�R5C�vpC�mz    ?z�HC� B���C��.C���    ?z�HC� B��C��WC��    ?z�HC� B��C�WC�    ?z�HC�� B�R�C�F[C�7�    ?z�HC� B��C�u�C�e�    ?z�HC� B�ļC��UC���    ?z�HC� B��`C���C���    ?z�HC� B�-�C���C��[    ?z�HC�� B�^�C� �C�%    ?z�HC�� B��C�F�C�/�    ?z�HC�� B�KC�j�C�S    ?z�HC�� B��C���C�t    ?z�HC�� B��C��4C���    ?z�HC�� B�1`C��C��D    ?z�HC�� B�R~C��bC��x    ?z�HC�� B�q]C��C��    ?z�HC�� B�C�IC���    ?z�HC�� B㨣C�6�C��    ?z�HC�� B��C�MC�'�    ?z�HD @ B��=C�a�C�:    ?z�HD � B��OC�t!C�J�    ?z�HD@ B��RC��>C�Z.    ?z�HD� B�0C��cC�g�    ?z�HD@ B��C���C�sf    ?z�HD� B�&�C���C�}�    ?z�HD@ B�0'C��.C���    ?z�HD� B�7�C���C��
    ?z�HD@ B�=<C��dC��a    ?z�HD� B�@�C���C���    ?z�HD@ B�BC��_C��    ?z�HD� B�AoC��UC��.    ?z�HD@ B�>�C�һC���    ?z�HD� B�:HC�ьC��    ?z�HD@ B�3�C���C��    ?z�HD� B�+HC��yC���    ?z�HD@ B� �C��{C��t    ?z�HD� B�nC���C�y�    ?z�HD	@ B�3C���C�n�    ?z�HD	� B���C��HC�b�    ?z�HD
@ B��C��=C�T�    ?z�HD
� B�ʃC���C�D�    ?z�HD@ B㲆C���C�3�    ?z�HD� B㘤C�o�C�!    ?z�HD@ B�}C�]�C��    ?z�HD� B�_nC�JC��    ?z�HD@ B�@0C�5 C��1    ?z�HD� B�C�bC���    ?z�HD@ B��C��C���    ?z�HD� B��jC��7C���    ?z�HD@ B��C��xC�u�    ?z�HD� B∻C��&C�W�    ?z�HD@ B�^�C���C�8�    ?z�HD� B�3&C�y�C��    ?z�HD@ B��C�YRC���    ?z�HD� B���C�7C�ҷ    ?z�HD@ B��C�}C��8    ?z�HD� B�s�C���C��x    ?z�HD@ B�?�C�ʉC�_�    ?z�HD� B�
 C��7C�6�    ?z�HD@ B���C�{C�3    ?z�HD� B���C�QdC���    ?z�HD@ B�__C�&�C��Z    ?z�HD� B�#dC��C���    ?z�HD@ B���C���C�X:    ?z�HD� Bߣ�C��C�':    ?z�HD@ B�azC�nC���    ?z�HD� B��C�='C���    ?z�HD@ B��C�
�C��G    ?z�HD� Bސ�C��C�W�    ?z�HD@ B�H�C��KC�!    ?z�HD� B���C�lgC��U    ?z�HD@ Bݳ?C�5VC���    ?z�HD� B�f~C��C�u�    ?z�HD@ B�?C�êC�:1    ?z�HD� B��zC��C���    ?z�HD@ B�woC�MEC��J    ?z�HD� B�$�C��C��a    ?z�HD@ B��3C���C�AE    ?z�HD� B�|C���C� �    ?z�HD@ B�%pC�S�C���    ?z�HD� B�ͫC��C�{�    ?z�HD@ B�t�C��"C�7�    ?z�HD� B��C��?C���    ?z�HD @ Bپ@C�H�C���    ?z�HD � B�aKC��C�e�    ?z�HD!@ B�C��dC�l    ?z�HD!� Bأ}C�u�C�Խ    ?z�HD"@ B�B�C�-�C���    ?z�HD"� B���C��sC�>�    ?z�HD#@ B�z�C��:C��,    ?z�HD#� B�!C�N�C���    ?z�HD$@ B֬rC��C�VY    ?z�HD$� B�C�C���C��    ?z�HD%@ B�ٟC�g�C���    ?z�HD%� B�n�C�C�eZ    ?z�HD&@ B�(C��.C�1    ?z�HD&� BԔ�C�x�C���    ?z�HD'@ B�&;C�&�C�l�    ?z�HD'� BӶ�C��`C��    ?z�HD(@ B�E�C��$C���    ?z�HD(� B��C�,�C�k�    ?z�HD)@ B�a-C��C�8    ?z�HD)� B��CC���C���    ?z�HD*@ B�xUC�+TC�c�    ?z�HD*� B�XC��C�	�    ?z�HD+@ BЋtC�{�C��7    ?z�HD+� B�fC�"�C�T,    ?z�HD,@ BϚ�C���C��    ?z�HD,� B� �C�nDC��)    ?z�HD-@ BΥ�C�C�=�    ?z�HD-� B�*	C��%C��`    ?z�HD.@ BͭgC�ZDC��t    ?z�HD.� B�/�C���C� �    ?z�HD/@ B̰�C���C���    ?z�HD/� B�/�C�?�C�^u    ?z�HD0@ BˮC���C��P    ?z�HD0� B�+�C��C���    ?z�HD1@ Bʨ!C��C�5�    ?z�HD1� B�#�C��C�т    ?z�HD2@ Bɞ�C�Z�C�l�    ?z�HD2� B��C���C�W    ?z�HD3@ Bȑ�C��_C��    ?z�HD3� B�
OC�02C�:n    ?z�HD4@ BǂC�ˁC��    ?z�HD4� B���C�e�C�j�    ?z�HD5@ B�n�C�  C�    ?z�HD5� B��*C��|C���    ?z�HD6@ B�X�C�2OC�/�    ?z�HD6� B�̏C���C��1    ?z�HD7@ B�?�C�b@C�Z    ?z�HD7� Bò"C���C���    ?z�HD8@ B�#�C��8C���    ?z�HD8� B�C�&NC�<    ?z�HD9@ B�/C���C��g    ?z�HD9� B�t�C�QC�<    ?z�HD:@ B��C��C���    ?z�HD:� B�RhC�y�C�_%    ?z�HD;@ B��C��C���    ?z�HD;� B�+C���C��r    ?z�HD<@ B���C�3KC��    ?z�HD<� B��C�ōC���    ?z�HD=@ B�`�C�W�C�*\    ?z�HD=� B��MC���C���    ?z�HD>@ B�+yC�y�C�D�    ?z�HD>� B���C�
BC�Ј    ?z�HD?@ B��C��ZC�\g    ?z�HD?� B�V�C�*C��    ?z�HD@@ B��7C��^C�r�    ?z�HD@� B�C�HC���    ?z�HDA@ B�|AC�֜C��;    ?z�HDA� B���C�d�C�o    ?z�HDB@ B�;�C��C��P    ?z�HDB� B��FC��5C�!�    ?z�HDC@ B���C�YC���    ?z�HDC� B�U%C��C�1q    ?z�HDD@ B���C�&�C���    ?z�HDD� B�C���C�?�    ?z�HDE@ B�j�C�>�C��d    ?z�HDE� B�ȑC��C�Ml    ?z�HDF@ B�&qC�U^C��L    ?z�HDF� B�� C��SC�ZB    ?z�HDG@ B���C�j�C���    ?z�HDG� B�/C��\C�b{    ?z�HDH@ B��hC�yC��    ?z�HDH� B���C�	SC�h:    ?z�HDI@ B�*�C���C���    ?z�HDI� B�~GC�MC�m    ?z�HDJ@ B�ѴC��cC��    ?z�HDJ� B�$�C�.XC�p�    ?z�HDK@ B�v�C���C��    ?z�HDK� B�ȺC�?fC�s�    ?z�HDL@ B�C���C��d    ?z�HDL� B�j�C�O�C�t�    ?z�HDM@ B��%C���C��i    ?z�HDM� B�
�C�_�C�uQ    ?z�HDN@ B�Z"C���C���    ?z�HDN� B��(C�n5C�s    ?z�HDO@ B��C���C��    ?z�HDO� B�;�C�|�C�n�    ?z�HDP@ B��[C�fC��D    ?z�HDP� B��C�� C�i�    ?z�HDQ@ B��C��C��	    ?z�HDQ� B�gC��3C�dm    ?z�HDR@ B���C��C��}    ?z�HDR� B�gC���C�a�    ?z�HDS@ B�U,C�*C��>    ?z�HDS� B��C��C�\~    ?z�HDT@ B���C�5�C���    ?z�HDT� B��C���C�N�    ?z�HDU@ B�L4C�AAC��    ?z�HDU� B�~NC���C~iY    ?z�HDV@ B��FC�L;C}Rm    ?z�HDV� B��C�ѥC|<4    ?z�HDW@ B��C�V�C{$�    ?z�HDW� B�V�C��/Cz�    ?z�HDX@ B��{C�a�Cx��    ?z�HDX� B��]C���Cw�|    ?z�HDY@ B��IC�k�Cv�N    ?z�HDY� B�4�C�8Cu��    ?z�HDZ@ B�jZC~�vCt��    ?z�HDZ� B���C}��Cs�     ?z�HD[@ B��C} �Cro7    ?z�HD[� B��C|CqX    ?z�HD\@ B�KmC{=Cp@�    ?z�HD\� B��&CzqCo)�    ?z�HD]@ B���Cy)�Cn#    ?z�HD]� B��Cx3�Cl��    ?z�HD^@ B�)�Cw>Ck�    ?z�HD^� B�\�CvHLCj��    ?z�HD_@ B��lCuR�Ci�[    ?z�HD_� B�ßCt\�Ch��    ?z�HD`@ B���CsgdCgk�    ?z�HD`� B�)'Crq�CfN+    ?z�HDa@ B�d�Cq|�Ce7)    ?z�HDa� B��Cp��Cd'l    ?z�HDb@ B��VCo��CcG    ?z�HDb� B�7�Cn�Cb
    ?z�HDc@ B���Cm��C`��    ?z�HDc� B��jCl�%C_�T    ?z�HDd@ B�dCk��C^��    ?z�HDd� B�\�Cj̿C]�o    ?z�HDe@ B���CiخC\�    ?z�HDe� B��RCh��C[�5    ?z�HDf@ B�7Cg�CZ��    ?z�HDf� B�w%Cf�CY�    ?z�HDg@ B���Cf
ACXj�    ?z�HDg� B��Ce�CWkN    ?z�HDh@ B�F�Cd%CVf    ?z�HDh� B�Cc2�CUV�    ?z�HDi@ B}��Cb@&CT>�    ?z�HDi� B|�CaM�CS&8    ?z�HDj@ Bz�lC`\CR9    ?z�HDj� By�C_j�CP��    ?z�HDk@ Bw�C^yiCO�    ?z�HDk� Bv5OC]�lCN�    ?z�HDl@ Bt��C\��CMѠ    ?z�HDl� BsPpC[�>CL�    ?z�HDm@ Bq�CZ�NCK��    ?z�HDm� Bp��CY��CJ��    ?z�HDn@ Bo��CX�gCI͇    ?z�HDn� Bnv�CW�	CH�    ?z�HDo@ Bmb�CV�CG��    ?z�HDo� BlLaCVCG2    ?z�HDp@ Bk3�CU"DCF.e    ?z�HDp� BjECT5 CEC�    ?z�HDq@ Bh�CSH�CDY�    ?z�HDq� Bg�BCR\UCCm�    ?z�HDr@ Bf��CQpsCB�x    ?z�HDr� Be�CP��CA��    ?z�HDs@ Bdv>CO��C@�[    ?z�HDs� BcOCN��C?��    ?z�HDt@ Bb&�CM�aC>ɴ    ?z�HDt� B`��CL�7C=�~    ?z�HDu@ B_��CK�C<�g    ?z�HDu� B^�
CKSC;��    ?z�HDv@ B]o3CJxC;�    ?z�HDv� B\;�CI6C:    ?z�HDw@ B[�CHM�C9V    ?z�HDw� BYτCGfBC8(�    ?z�HDx@ BX�QCF^C72�    ?z�HDx� BWZ�CE�{C6<�    ?z�HDy@ BV�CD��C5D�    ?z�HDy� BTܪCC�EC4Lq    ?z�HDz@ BS��CB��C3i>    ?z�HDz� BR�CB�C2��    ?z�HD{@ BQ�fCA5C1̯    ?z�HD{� BQ�C@< C0�    ?z�HD|@ BP.@C?Y�C0/    ?z�HD|� BOD�C>wdC/`    ?z�HD}@ BNY�C=��C.�    ?z�HD}� BMnFC<�EC-�2    ?z�HD~@ BL�wC;��C,�    ?z�HD~� BK�"C:�eC,!�    ?z�HD@ BJ�mC:�C+Q    ?z�HD� BI�lC944C*��    ?z�HD�  BH��C8U�C)��    ?z�HD�` BGԳC7w�C(�1    ?z�HD�� BF�C6�?C(x    ?z�HD�� BE��C5��C'>�    ?z�HD�  BD�xC4��C&nk    ?z�HD�` BD�C4LC%��    ?z�HD�� BC�C3(�C$��    ?z�HD�� BB*C2M�C#��    ?z�HD�  BA&C1slC#(�    ?z�HD�` B@/_C0��C"X>    ?z�HD�� B?7xC/�	C!��    ?z�HD�� B>>�C.�C ��    ?z�HD�  B=E�C.SC�    ?z�HD�` B<`�C-7�Cn    ?z�HD�� B;��C,aOC^R    ?z�HD�� B:��C+�
C�w    ?z�HD�  B9�?C*��C��    ?z�HD�` B9�C)��C'M    ?z�HD�� B8M�C)xCl!    ?z�HD�� B7}�C(8�C�~    ?z�HD�  B6�TC'e�C�N    ?z�HD�` B5�TC&�#C:�    ?z�HD�� B5�C%�9C�B    ?z�HD�� B4@�C$��Cź    ?z�HD�  B3o�C$1C    ?z�HD�` B2��C#OCP1    ?z�HD�� B1�C"�C��    ?z�HD�� B0�OC!��C��    ?z�HD�  B06�C �C'�    ?z�HD�` B/f�C �Co    ?z�HD�� B.�xCG�C��    ?z�HD�� B-ΘC{�C�    ?z�HD�  B,��C��CH�    ?z�HD�` B,,�C�C�    ?z�HD�� B+]�C	C�&    ?z�HD�� B*��CP+C$    ?z�HD�  B)�C��Co�    ?z�HD�` B(�C�SC�    ?z�HD�� B(%WC�YC�    ?z�HD�� B'lxC/7CX�    ?z�HD�  B&��ChQC�{    ?z�HD�` B%��C�IC
�3    ?z�HD�� B$�C��C
5    ?z�HD�� B$4C3C	��    ?z�HD�  B#L�CTC�I    ?z�HD�` B"y|C��C�    ?z�HD�� B!�{C��Cg&    ?z�HD�� B ��C�C�3    ?z�HD�  B pCJ.Cd    ?z�HD�` BM�C�fCY�    ?z�HD�� B�kC�HC�U    ?z�HD�� B�JC	�C J    ?z�HD�  B�CKCR�    ?z�HD�` B6�C��C�_    ?z�HD�� Bp�C�gC�s    ?z�HD�� B�|C�CP�    ?z�HD�  B�CVrC �>    ?z�HD�` B&�C
��B���    ?z�HD�� Bf~C	�1B���    ?z�HD�� B�UC	&B�[�    ?z�HD�  B�Cl�B��    ?z�HD�` B1�C��B�Ɗ    ?z�HD�� BuHC��B�|&    ?z�HD�� B��CDgB�2�    ?z�HD�  B�C��B��    ?z�HD�` BB�CטB���    ?z�HD�� B�C"=B�c_    ?z�HD�� B��Cm�B�+�    ?z�HD�  BXDC�B��    ?z�HD�` BT�C�B���    ?z�HD�� B��CSIB�M>    ?z�HD�� B�hC �EB�y    ?z�HD�  B,�B���B��T    ?z�HD�` Bs�B�~tB�a    ?z�HD�� B�vB�eB�N    ?z�HD�� B�5B���B�`    ?z�HD�  BB�B�bFB��}    ?z�HD�` B
��B��B�    ?z�HD�� B	�GB���B�V�    ?z�HD�� B	!XB�RjB��    ?z�HD�  Bm�B���B��    ?z�HD�` B�B�oB�    ?z�HD�� BcB�OB�~�    ?z�HD�� BV]B��,B�K�    ?z�HD�  B��B��B��    ?z�HD�` B��B�WjB���    ?z�HD�� BC-B��B۷o    ?z�HD�� B�{B�+Bڈ�    ?z�HD�  B�B�k�B�\7    ?z�HD�` B:"B� mB�2�    ?z�HD�� B��B��B�!    ?z�HD�� B �pB��B��    ?z�HD�  B ;QB�E�BԾ#    ?z�HD�` A�$�B��wBӘ�    ?z�HD�� A��QB�ZB�u�    ?z�HD�� A���B�w>B�SX    ?z�HD�  A�3B�56B�0�    ?z�HD�` A��B��*B��    ?z�HD�� A���Bݵ#B��|    ?z�HD�� A�K)B�v�B���    ?z�HD�  A���B�:JB˶?    ?z�HD�` A���B��Bʙ�    ?z�HD�� A�l�B���B�~�    ?z�HD�� A�$VB׋B�c�    ?z�HD�  A�݆B�TSB�KG    ?z�HD�` A�fB��B�3�    ?z�HD�� A�RuB��7B��    ?z�HD�� A��Bҵ�B��    ?z�HD�  A��bBу�B��    ?z�HD�` A�kB�S(B��h    ?z�HD�� A�I�B�#�B��s    ?z�HD�� A�
�B��lB���    ?z�HD�  A��=B��B���    ?z�HD�` A��B˝/B��    ?z�HD�� A�R
B�s	B���    ?z�HD�� A��B�I�B�u    ?z�HD�  A��KB�"NB�rc    ?z�HD�` A�BB��sB�h+    ?z�HD�� A�i9B��bB�^M    ?z�HD�� A�1�Bĳ�B�U�    ?z�HD�  A���BÑ�B�M�    ?z�HD�` A���B�p,B�G;    ?z�HD�� Aړ%B�P?B�B,    ?z�HD�� A�`�B�2�B�=�    ?z�HD�  A�/(B�B�:�    ?z�HD�` A���B���B�7�    ?z�HD�� A�ϚB�޶B�6�    ?z�HD�� AԢ	B���B�6�    ?z�HD�  A�u�B���B�8&    ?z�HD�` A�J-B��SB�:6    ?z�HD�� A��B��kB�>K    ?z�HD�� A���B�m�B�A�    ?z�HD�  A�βB�Z�B�G�    ?z�HD�` Aͧ�B�I�B�N    ?z�HD�� ÂRB�9�B�U�    ?z�HD�� A�^&B�*dB�^�    ?z�HD�  A�:�B��B�g�    ?z�HD�` A��B�MB�s�    ?z�HD�� A��aB�B��    ?z�HD�� A��B���B���    ?z�HD�  AŻB��B��S    ?z�HD�` AĞ2B��B��T    ?z�HD�� AÂYB��]B��    ?z�HD�� A�hB��B�ͺ    ?z�HD�  A�N�B��B��    ?z�HD�` A�6hB��NB��U    ?z�HD�� A��B��B�	�    ?z�HD�� A�
2B��dB�9    ?z�HD�  A���B��KB�6J    ?z�HD�` A��B��2B�O    ?z�HD�� A���B��vB�g�    ?z�HD�� A��KB��B��&    ?z�HD�  A��DB��B��    ?z�HD�` A���B���B��s    ?z�HD�� A���B���B��D    ?z�HD�� A��)B��B��^    ?z�HD�  A�|nB��B��    ?z�HD�` A�rjB��B�6�    ?z�HD�� A�iXB�+�B�W    ?z�HD�� A�bDB�:-B�|    ?z�HD�  A�[jB�J=B���    ?z�HD�` A�V�B�[�B�Û    ?z�HD�� A�R�B�m�B��.    ?z�HD�� A�OB���B�{    ?z�HD�  A�L�B���B�7�    ?z�HD�` A�L�B���B�`�    ?z�HD�� A�L�B��B���    ?z�HD�� A�NEB���B��v    ?z�HD�  A�Q_B��rB��&    ?z�HD�` A�U%B��B�    ?z�HD�� A�Z�B�*�B�<�    ?z�HD�� A�a
B�G@B�l    ?z�HD�  A�h�B�eyB���    ?z�HD�` A�p�B��{B��
    ?z�HD�� A�zuB���B��q    ?z�HD�� A�� B��pB�2,    ?z�HD�  A���B��bB�f<    ?z�HD�` A���B�SB��    ?z�HD�� A���B�1�B�р    ?z�HD�� A��#B�W�B�i    ?z�HD�  A��B�~�B�?�    ?z�HD�` A��aB��B~��    ?z�HD�� A��;B��JB}f�    ?z�HD�� A�;B��B{܊    ?z�HD�  A�(B�&�BzU�    ?z�HD�` A�0�B�S�Bx�|    ?z�HD�� A�G�B��	BwJ�    ?z�HD�� A�`B���Bu�a    ?z�HD�  A�{�B�ߋBtK    ?z�HD�` A���B�3Br�E    ?z�HD�� A��kB�CBqV�    ?z�HD�� A���B~�Bo��    ?z�HD�  A��B}U�Bnk    ?z�HD�` A�0B{�DBl�j    ?z�HD�� A�?�Bz.zBk��    ?z�HD�� A�c�Bx�FBj;    ?z�HD�  A���Bw8Bh��    ?z�HD�` A�� Bu��Bg=�    ?z�HD�� A��(Bs��Be��    ?z�HD�� A���Bro�Bdk�    ?z�HD�  A�'�Bp��Bcm    ?z�HD�` A�R"Boe?Ba�>    ?z�HD�� A�}�Bm�[B`=A    ?z�HD�� A���Blc�B^�.    ?z�HD�  A���Bj�B]|�    ?z�HD�` A�xBijUB\    ?z�HD�� A�4�Bg��BZÖ    ?z�HD�� A�eBfx�BYh�    ?z�HD�  A��RBeGBXP    ?z�HD�` A���Bc��BV��    ?z�HD�� A���BbmBUeQ    ?z�HD�� A�0WB`��BT�    ?z�HD�  A�e�B_@BR�:    ?z�HD�` A��B]ԛBQr.    ?z�HD�� A�_B\k�BP"�    ?z�HD�� A~"B[9BN�o    ?z�HD�  A|� BY�BM��    ?z�HD�` Az�yBX:~BLED    ?z�HD�� AyuBV�}BJ��    ?z�HD�� Aw�BUx�BI��    ?z�HD�  Avg6BT�BHs�    ?z�HD�` At�tBR�BG2C    ?z�HD�� Asb3BQd�BE��    ?z�HD�� Aq��BP�BD�@    ?z�HD�  ApdBN�OBCw�    ?z�HD�` An�lBMcBB<N    ?z�HD�� Amm�BL�BA�    ?z�HD�� Ak��BJ�!B?�e    ?z�HD�  Aj>BIq�B>��    ?z�HD�` Ai�BH%B=aY    ?z�HD�� Ag��BF�CB<-�    ?z�HD�� Af'�BE��B:��    ?z�HD�  Ad��BDJ�B9�g    ?z�HD�` AcLfBC�B8�Q    ?z�HD�� Aa�ABA�PB7r�    ?z�HD�� A`yfB@�B6H/    ?z�HD�  A_�B?A,B5    ?z�HD�` A]��B>�B3��    ?z�HD�� A\SMB<�B2֬    ?z�HD�� A[BB;�;B1��    ?z�HD�  AY�lB:T!B0��    ?z�HD�` AXkB9eB/��    ?z�HD�� AW �B7�B.r    ?z�HD�� AU�dB6�6B-[�    ?z�HD�  AT��B5��B,HH    ?z�HD�` ASH�B4USB+4W    ?z�HD�� ARHB3(B*!�    ?z�HD�� AP�4B1��B)�    ?z�HD�  AO{�B0�^B(�    ?z�HD�` AN:B/��B&�5    ?z�HD�� AL��B.��B%�    ?z�HD�� AK�B-]$B$ڌ    ?z�HD�  AJ}5B,9IB#�`    ?z�HD�` AIAuB+ B"ǀ    ?z�HD�� AHmB)��B!�{    ?z�HD�� AF�B(كB ��    ?z�HD�  AE��B'��B�    ?z�HD�` AD]�B&��B�    ?z�HD�� AC(�B%��B�,    ?z�HD�� AA�/B$o�B�    ?z�HD�  A@�+B#Z&B��    ?z�HD�` A?��B"E>B�M    ?z�HD�� A>_oB!1�B��    ?z�HD�� A=B�B  �B�    ?z�HD�  A<7=B�B��    ?z�HD�` A;+&B�B�    ?z�HD  A:�B�oB �    ?z�HD�� A9yB�Bl    ?z�HD�  A8^B�nB.V    ?z�HD�` A6�B�BE�    ?z�HDà A5�B�$B]&    ?z�HD�� A4�fB��Bu    ?z�HD�  A3�bB��B�c    ?z�HD�` A2�|BͱB�o    ?z�HDĠ A1�B�=B��    ?z�HD�� A0�B�$B܏    ?z�HD�  A/��B��B�{    ?z�HD�` A.~�B�B-    ?z�HDŠ A-qB��B0l    ?z�HD�� A,coB��B
MP    ?z�HD�  A+V�B��B	k�    ?z�HD�` A*J/B��B�    ?z�HDƠ A)>XB�B�I    ?z�HD�� A(36B�B�'    ?z�HD�  A'(�B �B��    ?z�HD�` A&�B
0�B�    ?z�HDǠ A%zB	CDB3�    ?z�HD�� A$=BVzBX�    ?z�HD�  A#UBjyB}M    ?z�HD�` A!��B��B��    ?z�HDȠ A �B��B Ǧ    ?z�HD�� AނB�A��!    ?z�HD�  AԖB��A�*�    ?z�HD�` A˒B�A�|A    ?z�HDɠ AâBHA���    ?z�HD�� A��B"�A�%�    ?z�HD�  A�-B BNA�}(    ?z�HD�` A��A��A���    ?z�HDʠ A�BA�<A�3�    ?z�HD�� A��A�T�A�    ?z�HD�  A�A��OA��	    ?z�HD�` A��A��3A�R�    ?z�HDˠ A�NA�<cA���    ?z�HD�� A�BA��A��    ?z�HD�  A��A���A�ix    ?z�HD�` AzCA�:CA���    ?z�HD̠ AhJA�|A��    ?z�HD�� AW�A��]A�v�    ?z�HD�  AH�A�P5A��s    ?z�HD�` A;KA�lA�6L    ?z�HD͠ A/iA�LA��    ?z�HD�� A%VA�}�A� �    ?z�HD�  A
}A��jA�jI    ?z�HD�` A	gA�R�A��o    ?z�HDΠ A�A���A�B�    ?z�HD�� A�A�2�Aز�    ?z�HD�  A	
AߦPA�#�    ?z�HD�` A�A�1Aՙ�    ?z�HDϠ AVAܔmA�`    ?z�HD�� A
6A�AҊ�    ?z�HD�  A�Aٌ�A��    ?z�HD�` A}A�Aτ,    ?z�HDР A �A֎�A��    ?z�HD�� @�I+A��Ả�    ?z�HD�  @�j�Aә�A�    ?z�HD�` @���A�"�Aɜ�    ?z�HDѠ @��!AЯ4A�*=    ?z�HD�� @���A�<�AƷ�    ?z�HD�  @��A��NA�H�    ?z�HD�` @�3(A�`HA��4    ?z�HDҠ @�`�A���A�n`    ?z�HD�� @AɌmA��    ?z�HD�  @��A�%�A���    ?z�HD�` @��hA��=A�7.    ?z�HDӠ @�-qA�`<A���    ?z�HD�� @�eA���A�o    ?z�HD�  @��A¢�A��    ?z�HD�` @���A�G�A��    ?z�HDԠ @��A��A�Q�    ?z�HD�� @�U�A��
A���    ?z�HD�  @ߕ�A�B�A��9    ?z�HD�` @���A��ZA�B�    ?z�HDՠ @�;A��pA��*    ?z�HD�� @�`�A�QQA���    ?z�HD�  @ئA�A�B    ?z�HD�` @��A���A��1    ?z�HD֠ @�7A�tvA���    ?z�HD�� @ӁpA�-�A�N    ?z�HD�  @��8A��A� ]    ?z�HD�` @��A���A���    ?z�HDנ @�gA�h�A�f    ?z�HD�� @̵�A�+�A��    ?z�HD�  @��A��?A���    ?z�HD�` @�U�A��'A��*    ?z�HDؠ @ǦpA��A�A�    ?z�HD�� @���A�J.A��D    ?z�HD�  @�I�A�A���    ?z�HD�` @�A���A�p�    ?z�HD٠ @��{A��3A�,2    ?z�HD�� @�B:A��A��    ?z�HD�  @���A�\.A��3    ?z�HD�` @���A�2A�d&    ?z�HDڠ @�;;A�
=A�"/    ?z�HD�� @��A��A���    ?z�HD�  @��A��*A��    ?z�HD�` @�1�A���A�^$    ?z�HD۠ @���A�}�A�    ?z�HD�� @���A�^�A��d    ?z�HD�  @� nA�BA��t    ?z�HD�` @�m�A�'DA�Y�    ?z�HDܠ @��A��A��    ?z�HD�� @�BA���A��V    ?z�HD�  @�JgA��A���    ?z�HD�` @���A�ͱA�L    ?z�HDݠ @�ѵA��$A��    ?z�HD�� @�
A���A���    ?z�HD�  @�MA��A�w!    ?z�HD�` @��+A��A�-_    ?z�HDޠ @��$A���A��N    ?z�HD�� @���A�{�A��S    ?z�HD�  @��A�s�A�?    ?z�HD�` @�:jA�lmA�Y    ?z�HDߠ @�Y�A�hfA} �    ?z�HD�� @�rgA�e8Azh�    ?z�HD�  @��qA�c�Aw�    ?z�HD�` @���A�dNAtؒ    ?z�HD� @���A�f6Ar �    ?z�HD�� @��:A�h�Ao_    ?z�HD�  @���A�o�Al�    ?z�HD�` @���A�x�Ak�    ?z�HD� @��-A���Ai]�    ?z�HD�� @���A���Ag��    ?z�HD�  @��NA���Ae��    ?z�HD�` @�A���AdP�    ?z�HD� @�	
A��-Ab��    ?z�HD�� @�[A�<Aa"    ?z�HD�  @�*A}�UA__D    ?z�HD�` @�&�A{�A]�,    ?z�HD� @�3�Az,�A\ �    ?z�HD�� @�A�Ax^WAZ��    ?z�HD�  @�RAv�AX�    ?z�HD�` @~�mAt��AWX(    ?z�HD� @|�1As
?AUŷ    ?z�HD�� @{AqHEAT4'    ?z�HD�  @yFAo�AR�n    ?z�HD�` @wt�Am��AQ�    ?z�HD� @u��Al[AO�1    ?z�HD�� @sݹAjb�AN>    ?z�HD�  @r,Ah��AL�    ?z�HD�` @pS�AgAK
p    ?z�HD� @n��AeU�AI��    ?z�HD�� @l�7Ac�@AH�    ?z�HD�  @kzAb4AF��    ?z�HD�` @i`�A`cAE!�    ?z�HD� @g��A^�.AC�    ?z�HD�� @e�/A]%�AB>/    ?z�HD�  @dI�A[�EA@��    ?z�HD�` @b�AY�>A?bb    ?z�HD� @`�AXZ�A=��    ?z�HD�� @_LAVɋA<��    ?z�HD�  @]��AU8VA;.L    ?z�HD�` @\�AS��A9��    ?z�HD� @Zh�AR"A8l�    ?z�HD�� @X�FAP�-A7<    ?z�HD�  @W3�AOfA5�_    ?z�HD�` @U�kAM�A4\*    ?z�HD� @T�AL�A3,    ?z�HD�� @RyxAJ��A1�    ?z�HD�  @P�AI�A0`�    ?z�HD�` @O_sAG��A/Y    ?z�HD� @M�AF+�A-�'    ?z�HD�� @LO�AD�uA,x�    ?z�HD�  @J̃ACIKA+0G    ?z�HD�` @IK>AAڴA)�    ?z�HD� @G̯A@oBA(��    ?z�HD�� @FQOA?�A'c�    ?z�HD�  @D׼A=�<A&$     ?z�HD�` @Ca�A<<�A$�    ?z�HD�� @BlA:�:A#��    ?z�HD�� @@��A9}^A"��    ?z�HD�  @?��A8!�A!�5    ?z�HD�` @>l]A6ǅA ��    ?z�HD� @=:A5p�At�    ?z�HD�� @<GA4�Af    ?z�HD�  @:��A2��AX�    ?z�HD�` @9��A1x�AM)    ?z�HD� @8y�A0+`AAG    ?z�HD�� @7M=A.�/A7�    ?z�HD�  @6!�A-��A.�    ?z�HD�` @4�1A,PA'p    ?z�HD� @3ΦA+3A".    ?z�HD�� @2�bA)��A0    ?z�HD�  @1��A(�|A�    ?z�HD�` @0Z�A'J�A8    ?z�HD� @/7A&mAG    ?z�HD�� @.3A$ԊAD    ?z�HD�  @,��A#�A�    ?z�HD�` @+ҹA"h�A    ?z�HD� @*��A!5A#�    ?z�HD�� @)�;A �A*8    ?z�HD�  @(z-A��A1�    ?z�HD�` @'_0A�oA;%    ?z�HD� @&FA}�AEL    ?z�H
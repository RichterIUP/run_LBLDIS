CDF     
      n_wnum         n_cloud_layers     
   n_instances       n_height   E         Author        iLBLDIS developed by D.D. Turner, Space Science and Engineering Center, UW-Madison, dturner@ssec.wisc.edu       Version       G$Id: lblrtm_disort.c,v 1.44 2014/01/13 21:23:24 merrelli Release_3_0 $     Number_of_streams         16     RT_solver_number      0      Calculation_zenith_angle      180.000000 degrees      Calculation_zenith_angle_comment      I0 degrees is nadir (upwelling) while 180 degrees is zenith (downwelling)       
Input_flag        1      Path_to_LBLRTM_gas_ODs        {/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/run_LBLDIS/out/.lblrtm_node52.cluster_20210813_16064141CEST     SSP_database_number_0         m/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/ssp_database/ssp_db.mie_wat.gamma_sigma_0p100       SSP_database_number_1         m/home/phi.richter/retrieval_recode/Total_Cloud_Water_retrieval/ssp_database/ssp_db.mie_ice.gamma_sigma_0p100       Wavenumber_comment        &Using selected microwindows option #0      Phase_function_comment        Using real phase functions     Solar_zenith_angle        2No solar source input included in the calculation      Solar_azimuth_angle       2No solar source input included in the calculation      Sun_earth_distance        2No solar source input included in the calculation      Solar_source_datafile         2No solar source input included in the calculation            wnum                	long_name         wavenumber     units         cm-1            X   radiance                   	long_name         Computed radiance      units         mW / (m2 sr cm-1)           \   flux_up                    	long_name         Computed upwelling flux    units         mW / (m2 cm-1)          `   	flux_down                      	long_name         "Computed downwelling diffuse flux      units         mW / (m2 cm-1)          d   	flux_beam                      	long_name         Computed direct beam flux      units         mW / (m2 cm-1)          h   
cld_height                 	long_name         Cloud height       units         km        (  L   	cld_dbnum                  	long_name         Cloud database number      units         	unitless       comment       0database number of SSP used for the cloud layer         t   cld_tau                   	long_name         Cloud optical depth    units         	unitless          (  �   cld_tau_wnum               	long_name         1Reference wavenumber for the cloud optical depth       units         cm-1       comment       ZA value of -1 implies the optical depth is the geometric limit value (i.e., where Qe = 2)         (  �   cld_reff               	long_name          Cloud particle effective radius    units         microns       (  �   sfc_temperature              	long_name         Surface temperature    units         K      comment       oA value of -1 indicates that the surface temperature is assumed to be the same as the lowest level temperature              sfc_emissivity                  	long_name         Surface emissivity     units         	unitless            l   pwv              	long_name         Precipitable water vapor       units         cm             height                 	long_name         Level height       units         km            pressure               	long_name         Level pressure     units         mb            temperature                	long_name         Level temperature      units         K          0   mixing_ratio               	long_name         Level water vapor mixing ratio     units         g/kg       comment       �This profile is estimated from the number of molecules per layer that is in the LBLRTM TAPE7 output.  This profile should be made to agree with the precipitable water vapor (PWV) amount that was also derived from the TAPE7, as the PWV is more accurate.           D=���=���>L��>L��>���>���>���>���?   ?                  @       @       @       @       @       ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  @�  A�  @�  A�  @�  A�  @�  A�  @�  A�  ��      <�t�=u=���>��>L��>�  >���>�33>�Q�>�p�>\>Ǯ>���>�ff?   ?�?�?�R?8Q�?G�?O\)?\(�?k�?���?�  ?�33?���?�ff@   @��@��@&ff@333@@  @L��@fff@s33@�  @�ff@���@�ff@�  @���@�33@���@�ff@�  @ٙ�@�33@���@�ffA   AffA��A33A��A   A(  A0  A8  A@  AH  AP  A`  Ap  A�  A�  A�  A�  D}hRD|{Dz��Dy@�Dw�HDv$Dt��Ds�Dr�Drr�Dr$�DqևDq��Dp�Dn��Dm�DlyyDk+Dh,JDfyHDe��Dd:�Db�D]�RDY��DU�HDPsdDKb=DFi�DA�yD<D8�D3��D/
-D*�bD"&VDD�PD9D.�D��D�D�9C�o�C�RC��CܤCӤZC��}CC�`bC���C�fFC��9C���C���C�G+Ct��Ca�CP�)CALC2�1C%��C[�B��FB�wLBi=qA��AQ��C���C���C�Y�C�&fC��3C�ٚC��3C�� C�s3C�s3C�s3C�s3C�s3C�&fC�� C�&fC�33C��3C�� C��3C��3C�� C��3C��fC���C�s3C��C�ffC��3C�  C�s3C��C���C�  C�Y�C�  C�Y�C���C�Y�C��fC��C�33C�33C~33C{�fCy��CwL�Cu�Cr�fCp��Cn�3ClffCi33Ce��Cb�C^��C[��CXL�CVL�CX��C\��C_� C`L�C`�fCaL�Cb� Cc33Cc�fCk�f                                                                                                                                                                                                                                                                                    |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  |�  C� B�	"C��%C�oi    ?z�HC� B�g�C��C��h    ?z�HC� B���C�7C��    ?z�HC� B�}C�_C�I3    ?z�HC� B�t�C��bC��    ?z�HC� B��VC�	nC��    ?z�HC� B��C�K�C�H    ?z�HC� B�kfC��C�P�    ?z�HC� B��C���C��    ?z�HC� B��C�!C��    ?z�HC� B�LJC�A�C�&    ?z�HC� B�vC�z_C�9�    ?z�HC� B��0C��oC�n�    ?z�HC� B��C��C���    ?z�HC� B�V�C�FC��    ?z�HC�� B�;C�JqC�-    ?z�HC� B�͋C�z	C�2G    ?z�HC� B�wC��DC�]�    ?z�HC� B�;C���C��}    ?z�HC� B�nSC���C���    ?z�HC�� B�NC�$�C��H    ?z�HC�� B���C�JC���    ?z�HC�� B��ZC�oC��    ?z�HC�� B�$nC���C�@�    ?z�HC�� B�LMC��UC�`�    ?z�HC�� B�q�C��gC�~L    ?z�HC�� B唷C��C���    ?z�HC�� B�C�	�C���    ?z�HC�� B���C�#8C��%    ?z�HC�� B��AC�;6C��    ?z�HC�� B�
jC�Q=C���    ?z�HD @ B�"_C�e�C�
�    ?z�HD � B�8C�xvC�<    ?z�HD@ B�K�C��`C�+z    ?z�HD� B�])C���C�9�    ?z�HD@ B�lxC��{C�E|    ?z�HD� B�y�C���C�QK    ?z�HD@ B愺C���C�Y�    ?z�HD� B��C��(C�a!    ?z�HD@ B攪C���C�f�    ?z�HD� B�zC���C�j�    ?z�HD@ B�AC�ԀC�mg    ?z�HD� B�C��wC�m�    ?z�HD@ B曵C��AC�m    ?z�HD� B�nC�ծC�j�    ?z�HD@ B�C���C�g    ?z�HD� B��C�ΙC�a^    ?z�HD@ B�zC�ȳC�Z`    ?z�HD� B�w%C��:C�Q�    ?z�HD	@ B�i�C��.C�G�    ?z�HD	� B�Z�C���C�<g    ?z�HD
@ B�ITC���C�.�    ?z�HD
� B�6C���C��    ?z�HD@ B� �C���C��    ?z�HD� B�	�C�t!C��m    ?z�HD@ B���C�b
C��    ?z�HD� B���C�N}C��Q    ?z�HD@ B�ZC�9�C���    ?z�HD� B��C�#C��    ?z�HD@ B�z�C�C���    ?z�HD� B�X�C��C�t�    ?z�HD@ B�4�C��0C�X�    ?z�HD� B�C���C�;M    ?z�HD@ B��C���C��    ?z�HD� B�yC�~LC��g    ?z�HD@ B䓐C�]�C��Q    ?z�HD� B�gC�<7C���    ?z�HD@ B�8�C�6C���    ?z�HD� B��C��C�oO    ?z�HD@ B��+C���C�H�    ?z�HD� B��C��QC� �    ?z�HD@ B�n�C��C��l    ?z�HD� B�8MC�VC���    ?z�HD@ B� %C�+<C���    ?z�HD� B��OC��/C�t!    ?z�HD@ B��C�юC�E$    ?z�HD� B�M�C���C��    ?z�HD@ B��C�s(C��    ?z�HD� B�ΝC�A�C���    ?z�HD@ B��C�?C��    ?z�HD� B�IpC��C�J�    ?z�HD@ B��C��hC�    ?z�HD� B�>C�q�C�ު    ?z�HD@ B�v�C�:C��K    ?z�HD� B�-BC�5C�m;    ?z�HD@ B��C��bC�2F    ?z�HD� BߖaC���C��F    ?z�HD@ B�H�C�RbC��    ?z�HD� B���C�jC�}    ?z�HD@ BީyC�שC�>�    ?z�HD� B�W�C���C���    ?z�HD@ B��C�X�C��    ?z�HD� BݰjC��C�|:    ?z�HD@ B�Z�C�աC�8�    ?z�HD� B��C��C���    ?z�HD @ Bܫ3C�N4C���    ?z�HD � B�Q�C��C�i�    ?z�HD!@ B���C�°C�"�    ?z�HD!� Bۚ�C�{C�ڔ    ?z�HD"@ B�=#C�3C���    ?z�HD"� B��SC��C�G~    ?z�HD#@ B�~0C���C��    ?z�HD#� B��C�T}C���    ?z�HD$@ Bٺ?C��C�c�    ?z�HD$� B�VvC��4C��    ?z�HD%@ B��C�m4C��    ?z�HD%� B؋~C�%C�w    ?z�HD&@ B�$LC��LC�&�    ?z�HD&� B׻�C�~C�Տ    ?z�HD'@ B�R�C�,�C��'    ?z�HD'� B��C��
C�/�    ?z�HD(@ B�|iC���C�۬    ?z�HD(� B��C�2�C���    ?z�HD)@ BբC�݃C�15    ?z�HD)� B�34C���C��V    ?z�HD*@ B��UC�0�C��9    ?z�HD*� B�R�C�كC�*�    ?z�HD+@ B���C���C��%    ?z�HD+� B�m�C�(�C�x�    ?z�HD,@ B��.C���C��    ?z�HD,� B҅bC�t]C��    ?z�HD-@ B��C�C�f�    ?z�HD-� Bљ"C��C�
L    ?z�HD.@ B�!�C�`?C���    ?z�HD.� BЩ0C�C�N�    ?z�HD/@ B�/�C��hC���    ?z�HD/� BϵC�E�C���    ?z�HD0@ B�9�C��>C�/�    ?z�HD0� Bν�C���C���    ?z�HD1@ B�@�C�%C�l�    ?z�HD1� B�C�ÑC�
�    ?z�HD2@ B�DC�a<C��,    ?z�HD2� B��uC��C�DR    ?z�HD3@ B�DC���C��     ?z�HD3� B���C�6JC�{i    ?z�HD4@ B�A'C�љC��    ?z�HD4� BʾqC�l�C���    ?z�HD5@ B�;C�?C�IF    ?z�HD5� Bɶ�C���C��    ?z�HD6@ B�2C�8�C�z|    ?z�HD6� BȬvC��C�%    ?z�HD7@ B�&1C�i=C���    ?z�HD7� Bǟ7C� 7C�@�    ?z�HD8@ B��C��3C��     ?z�HD8� BƏHC�--C�l     ?z�HD9@ B�QC���C� �    ?z�HD9� B�|�C�X2C���    ?z�HD:@ B��nC���C�*    ?z�HD:� B�g�C���C��R    ?z�HD;@ B��C��C�P8    ?z�HD;� B�O�C���C��'    ?z�HD<@ B��!C�:�C�t�    ?z�HD<� B�4C���C�}    ?z�HD=@ B��lC�^�C��<    ?z�HD=� B�*C��+C�'�    ?z�HD>@ B��4C��?C���    ?z�HD>� B���C��C�G�    ?z�HD?@ B�eC���C��1    ?z�HD?� B�ӏC�1VC�e�    ?z�HD@@ B�AzC���C��B    ?z�HD@� B���C�O�C���    ?z�HDA@ B��C�ެC�    ?z�HDA� B��_C�l�C���    ?z�HDB@ B��1C���C�*�    ?z�HDB� B�_gC��xC��C    ?z�HDC@ B��AC��C�C>    ?z�HDC� B�4�C��]C��4    ?z�HDD@ B��tC�.�C�Z�    ?z�HDD� B��C��,C��    ?z�HDE@ B�qC�GC�q    ?z�HDE� B���C�ҷC���    ?z�HDF@ B�B<C�^
C��    ?z�HDF� B��C��C�"    ?z�HDG@ B�PC�s�C���    ?z�HDG� B�w�C��VC�#+    ?z�HDH@ B���C���C��-    ?z�HDH� B�CwC��C�4�    ?z�HDI@ B���C��CC��E    ?z�HDI� B�vC�%�C�Er    ?z�HDJ@ B�q�C���C��    ?z�HDJ� B��C�8-C�U    ?z�HDK@ B�9�C���C��<    ?z�HDK� B��-C�IkC�c�    ?z�HDL@ B� =C���C��    ?z�HDL� B�b�C�Y�C�q�    ?z�HDM@ B��PC���C��    ?z�HDM� B�'RC�i�C�~L    ?z�HDN@ B��C��C�h    ?z�HDN� B��C�y8C���    ?z�HDO@ B�I�C� ]C��    ?z�HDO� B��C���C��'    ?z�HDP@ B�	�C��C�M    ?z�HDP� B�iRC��pC��K    ?z�HDQ@ B�ȑC�0C�$     ?z�HDQ� B�'�C���C���    ?z�HDR@ B���C�)QC�-r    ?z�HDR� B���C���C��	    ?z�HDS@ B�D�C�5�C�6v    ?z�HDS� B���C��*C���    ?z�HDT@ B���C�B>C�>S    ?z�HDT� B�[*C��:C���    ?z�HDU@ B��C�NC�C�    ?z�HDU� B�#C���C��Z    ?z�HDV@ B�fMC�Y�C�H�    ?z�HDV� B��TC��WC���    ?z�HDW@ B��C�d�C�M    ?z�HDW� B�p�C��SC��    ?z�HDX@ B���C�o�C~��    ?z�HDX� B�!C���C}��    ?z�HDY@ B�x�C�z�C|�7    ?z�HDY� B��C� `C{�N    ?z�HDZ@ B�'C�Cz�$    ?z�HDZ� B�}�C~DCy��    ?z�HD[@ B�աC} �Cx��    ?z�HD[� B�-,C|+�Cw��    ?z�HD\@ B��wC{6gCv��    ?z�HD\� B�ۺCzACu�b    ?z�HD]@ B�2 CyK�Ct��    ?z�HD]� B��_CxV�Cs��    ?z�HD^@ B���CwabCr��    ?z�HD^� B�4<Cvl7Cq�g    ?z�HD_@ B��iCuwCp��    ?z�HD_� B��NCt�Co�    ?z�HD`@ B�2�Cs�Cn�-    ?z�HD`� B��ICr�:Cm�    ?z�HDa@ B��cCq�|Cl��    ?z�HDa� B�8�Cp��Ck�N    ?z�HDb@ B���Co�kCj��    ?z�HDb� B��Cn�Ciښ    ?z�HDc@ B�G�Cm��Chߋ    ?z�HDc� B���Cl��Cg�    ?z�HDd@ B��Ck��Cf�    ?z�HDd� B�V�Cj�2Ce��    ?z�HDe@ B���Cj�Cd�    ?z�HDe� B�-Ci;Cc�    ?z�HDf@ B�d�ChCb�N    ?z�HDf� B���Cg)Cb     ?z�HDg@ B��Cf63C`��    ?z�HDg� B�nrCeC�C`�    ?z�HDh@ B���CdQvC_�    ?z�HDh� B�(�Cc_fC^r    ?z�HDi@ B��Cbm�C]�    ?z�HDi� B��kCa{�C\!�    ?z�HDj@ B�4�C`��C[&�    ?z�HDj� B��C_��CZ.    ?z�HDk@ B��C^�CY3�    ?z�HDk� B�GC]�HCX:�    ?z�HDl@ B��	C\��CWAT    ?z�HDl� B���C[�jCVGs    ?z�HDm@ B�WCZ�zCUM�    ?z�HDm� B��CY�TCT[u    ?z�HDn@ B�(wCY
�CSo�    ?z�HDn� B���CX�CR�b    ?z�HDo@ B��CW.=CQ��    ?z�HDo� B�oWCV@�CP�    ?z�HDp@ B���CUSCO��    ?z�HDp� B�G�CTfKCN�m    ?z�HDq@ B���CSy�CM�    ?z�HDq� B��CR��CL�9    ?z�HDr@ BfCQ��CL�    ?z�HDr� B}�CP�CK%�    ?z�HDs@ B|�1CO�yCJ9�    ?z�HDs� B{��CN��CINm    ?z�HDt@ Bzm�CM�lCHb�    ?z�HDt� ByB�CM�CGw    ?z�HDu@ BxCL""CF�d    ?z�HDu� Bv��CK8CE��    ?z�HDv@ Bu�SCJPkCD�A    ?z�HDv� Bt�:CIg�CCɘ    ?z�HDw@ Bsc�CH�5CB�    ?z�HDw� Br5�CG�6CA�    ?z�HDx@ Bq�CF��CAL    ?z�HDx� Bo��CE�8C@�    ?z�HDy@ Bn�^CD�C?0<    ?z�HDy� Bmw�CC��C>D�    ?z�HDz@ BlYDCCC=`Y    ?z�HDz� BkK�CB5�C<��    ?z�HD{@ Bj>�CAQ�C;�    ?z�HD{� Bi0�C@nxC:�L    ?z�HD|@ Bh#�C?�C9��    ?z�HD|� BgC>�WC9.    ?z�HD}@ Bf+C=ǒC869    ?z�HD}� Bd��C<��C7[     ?z�HD~@ Bc��C<6C6Y    ?z�HD~� Bb߷C;$�C5�    ?z�HD@ Ba�5C:DfC4ʚ    ?z�HD� B`��C9d�C3�    ?z�HD�  B_�5C8��C3�    ?z�HD�` B^�&C7��C2=    ?z�HD�� B]��C6��C1c�    ?z�HD�� B\��C5��C0��    ?z�HD�  B[��C5�C/��    ?z�HD�` BZv(C43�C.��    ?z�HD�� BYifC3XC.S    ?z�HD�� BX])C2|�C-*X    ?z�HD�  BWP�C1��C,R    ?z�HD�` BVD�C0��C+{j    ?z�HD�� BU8�C/��C*��    ?z�HD�� BT,�C/xC)�i    ?z�HD�  BS!�C.=XC(�S    ?z�HD�` BR%C-e�C(&�    ?z�HD�� BQ"(C,�	C'W�    ?z�HD�� BP&�C+�C&�B    ?z�HD�  BO+*C*�LC%�Q    ?z�HD�` BN0�C*�C$�    ?z�HD�� BM6�C)9C$#,    ?z�HD�� BL=1C(eC#W    ?z�HD�  BKDC'��C"��    ?z�HD�` BJKeC&��C!��    ?z�HD�� BIS=C%�C �    ?z�HD�� BH[MC%�C +�    ?z�HD�  BGc�C$I�Cb>    ?z�HD�` BFl%C#y}C��    ?z�HD�� BEu<C"��C�6    ?z�HD�� BDaC!ڄC1    ?z�HD�  BC��C!�CAO    ?z�HD�` BB��C >Cz>    ?z�HD�� BA��Cp�C��    ?z�HD�� B@�TC�C�3    ?z�HD�  B?�YC�C(�    ?z�HD�` B>�C�Cc�    ?z�HD�� B=ծCA�C�    ?z�HD�� B<�Cw�C�U    ?z�HD�  B;�{C�*CM    ?z�HD�` B;nC�JCUm    ?z�HD�� B:�CC��    ?z�HD�� B9(�CU�C�w    ?z�HD�  B86C��CG    ?z�HD�` B7I"C�KCR�    ?z�HD�� B6]�C�C��    ?z�HD�� B5q�C=�CՈ    ?z�HD�  B4��Cy]CC    ?z�HD�` B3�AC��CY    ?z�HD�� B2��C�C��    ?z�HD�� B1ƤC0NC�d    ?z�HD�  B0��Cn�C%     ?z�HD�` B/�C��Cj�    ?z�HD�� B/�C�5C��    ?z�HD�� B.,�C-�C
�]    ?z�HD�  B-GUCnsC
>g    ?z�HD�` B,c^C�C	�P    ?z�HD�� B+�C�XCΫ    ?z�HD�� B*��C5KC�    ?z�HD�  B)��Cx�Cac    ?z�HD�` B(ٮC
�2C��    ?z�HD�� B'��C
(C��    ?z�HD�� B'^C	G�CB�    ?z�HD�  B&<HC�C�J    ?z�HD�` B%^�C�C܃    ?z�HD�� B$�C�C*;    ?z�HD�� B#�CeCx�    ?z�HD�  B"�hC�	C�2    ?z�HD�` B!�C��CO    ?z�HD�� B!KCBC iY    ?z�HD�� B C�C�B�wt    ?z�HD�  Bt�C��B�"8    ?z�HD�` B��C%*B���    ?z�HD�� B�mCr9B�b�    ?z�HD�� B��C ��B�
�    ?z�HD�  BC ^B���    ?z�HD�` B:�B���B�]Z    ?z�HD�� Bf�B�ZpB�6    ?z�HD�� B�B��QB��D    ?z�HD�  B��B���B�a6    ?z�HD�` B��B�@�B�'    ?z�HD�� B�B��1B���    ?z�HD�� BM�B���B�re    ?z�HD�  BB�38B�#R    ?z�HD�` B�eB��B��    ?z�HD�� B�2B�B��    ?z�HD�� B�B�3�B�F     ?z�HD�  BJ}B���B���    ?z�HD�` BB�FB綳    ?z�HD�� B��B�>�B�pA    ?z�HD�� B�)B��B�-�    ?z�HD�  B"sB��B��~    ?z�HD�` B[~B�U�B⪽    ?z�HD�� B�fB��B�j�    ?z�HD�� BϭB�°B�-�    ?z�HD�  B
�B�y�B��    ?z�HD�` BFxB�4+Bݵ�    ?z�HD�� B
�zB��|B�z�    ?z�HD�� B	��B�aB�B�    ?z�HD�  B��B�g�B�
z    ?z�HD�` B=�B�&�B�ӄ    ?z�HD�� B}B���Bמ�    ?z�HD�� B��BܨB�jT    ?z�HD�  B�B�jGB�81    ?z�HD�` B@�B�.B��    ?z�HD�� B�BB��\B�׶    ?z�HD�� BƳB׻�BѨO    ?z�HD�  B
�BքPB�{u    ?z�HD�` BPB�K8B�No    ?z�HD�� B�iB�iB�#�    ?z�HD�� B �SB��
B��<    ?z�HD�  B #(Bѱ<B�Ђ    ?z�HD�` A��PBЀBʪ�    ?z�HD�� A�gbB�PYBɂ�    ?z�HD�� A��kB�"CB�^�    ?z�HD�  A��oB���B�;K    ?z�HD�` A�$[B�ǫB��    ?z�HD�� A���BʞMB���    ?z�HD�� A�TB�s^B��&    ?z�HD�  A��sB�M&B¹�    ?z�HD�` A�B�&�B��    ?z�HD�� A�'B� �B��    ?z�HD�� A��B��B�f    ?z�HD�  A�eaBú�B�K;    ?z�HD�` A�DBB�3�    ?z�HD�� A�cB�y'B�|    ?z�HD�� A�O;B�Y:B��    ?z�HD�  A��tB�<tB��    ?z�HD�` A��B��B��B    ?z�HD�� A�F{B�B��s    ?z�HD�� A��fB��jB��p    ?z�HD�  A䝘B��\B���    ?z�HD�` A�L-B���B��{    ?z�HD�� A��6B��dB��L    ?z�HD�� A�vB��|B��`    ?z�HD�  A�^B�&B�zT    ?z�HD�` A��B�mbB�qL    ?z�HD�� A��[B�]1B�g�    ?z�HD�� A�}�B�M B�_�    ?z�HD�  A�6?B�A�B�Z�    ?z�HD�` A���B�4}B�U�    ?z�HD�� A׫�B�)B�R�    ?z�HD�� A�hB��B�Q_    ?z�HD�  A�&]B�9B�N�    ?z�HD�` A��B��B�Q2    ?z�HD�� Aҧ`B�'B�QI    ?z�HD�� A�jEB��B�U/    ?z�HD�  A�.�B��_B�Y$    ?z�HD�` A���B��XB�]    ?z�HD�� AͻJB��PB�cB    ?z�HD�� ĀHB���B�kq    ?z�HD�  A�NB���B�s�    ?z�HD�` A�;B��<B�}    ?z�HD�� A��PB� �B��@    ?z�HD�� AǶ1B�KB���    ?z�HD�  AƅmB�B��]    ?z�HD�` A�WZB�pB���    ?z�HD�� A�)�B�kB��V    ?z�HD�� A��SB�$gB��:    ?z�HD�  A��<B�.�B��    ?z�HD�` A��8B�;�B��    ?z�HD�� A���B�H�B�
Q    ?z�HD�� A�]�B�V�B��    ?z�HD�  A�9�B�f�B�8    ?z�HD�` A�6B�xB�N�    ?z�HD�� A��gB���B�j    ?z�HD�� A��eB��B��s    ?z�HD�  A��zB��B��I    ?z�HD�` A���B�ƝB��/    ?z�HD�� A�{�B���B���    ?z�HD�� A�`�B���B��a    ?z�HD�  A�GYB�2B��    ?z�HD�` A�/sB�(�B�7    ?z�HD�� A�B�D�B�X{    ?z�HD�� A�FB�a6B�|0    ?z�HD�  A��DB�~�B��<    ?z�HD�` A��B���B��~    ?z�HD�� A�́B���B��    ?z�HD�� A��<B��VB��    ?z�HD�  A���B��B�;�    ?z�HD�` A��B�%�B�f{    ?z�HD�� A���B�JB��Q    ?z�HD�� A��]B�pB��u    ?z�HD�  A���B���B��L    ?z�HD�` A�|pB���B�	    ?z�HD�� A�v'B���B�G�    ?z�HD�� A�q�B��B�w�    ?z�HD�  A�m�B�=�B���    ?z�HD�` A�k�B�j�B�ۃ    ?z�HD�� A�j�B��B��    ?z�HD�� A�k;B�ƙB~�U    ?z�HD�  A�nB��B|�    ?z�HD�` A�rWB�&qB{c    ?z�HD�� A�wVB�YBy�2    ?z�HD�� A�}�B�BxBZ    ?z�HD�  A���B}�Bv��    ?z�HD�` A��sB{�*Bu/�    ?z�HD�� A���BzX�Bs�i    ?z�HD�� A��Bx�/Br!a    ?z�HD�  A��pBw6�Bp�v    ?z�HD�` A���Bu��Bo&    ?z�HD�� A��!Bt(Bm��    ?z�HD�� A��|Br�dBl �    ?z�HD�  A��Bq3Bj��    ?z�HD�` A��Bo��Bi,    ?z�HD�� A�'Bn	�Bg�W    ?z�HD�� A�/�Bl��Bf?    ?z�HD�  A�FYBkPBd̫    ?z�HD�` A�^�Bi��Bc[�    ?z�HD�� A�w�Bh�Ba��    ?z�HD�� A��aBf��B`~T    ?z�HD�  A���Be&xB_�    ?z�HD�` A�˻Bc�(B]�    ?z�HD�� A���Bb?�B\B�    ?z�HD�� A�	ZB`ϮBZ��    ?z�HD�  A�)�B_b�BYx�    ?z�HD�` A�LB]��BXN    ?z�HD�� A�oIB\�XBV�    ?z�HD�� A��mB[%BUZ{    ?z�HD�  A���BY�SBS��    ?z�HD�` A��HBXZ�BR��    ?z�HD�� A��BV��BQLR    ?z�HD�� A�/�BU�PBO��    ?z�HD�  A�Y�BT9�BN�P    ?z�HD�` A��FBRݾBMO�    ?z�HD�� A���BQ��BK�    ?z�HD�� A�3BP+
BJ��    ?z�HD�  A~�BNՠBIc^    ?z�HD�` A||IBM�4BHr    ?z�HD�� Az��BL/�BFλ    ?z�HD�� AyB^BJ��BE�[    ?z�HD�  Aw��BI�BDBN    ?z�HD�` Av~BHB#BB�)    ?z�HD�� At|�BF��BA��    ?z�HD�� Ar��BE�TB@|�    ?z�HD�  AqZBDfB?=�    ?z�HD�` Ao�XBC!�B>    ?z�HD�� An?_BA��B<�M    ?z�HD�� Al�B@�PB;��    ?z�HD�  Ak.]B?[�B:W�    ?z�HD�` Ai��B>B9!J    ?z�HD�� Ah(6B<�IB7�    ?z�HD�� Af�4B;��B6�_    ?z�HD�  Ae0�B:nB5�    ?z�HD�` Ac�DB98B4d�    ?z�HD�� AbB�B81B38�    ?z�HD�� A`��B6��B2x    ?z�HD�  A_\	B5��B0�F    ?z�HD�` A]�B4l�B/��    ?z�HD�� A\}B3>�B.�v    ?z�HD�� A[rB2UB-|�    ?z�HD�  AY��B0�dB,Y�    ?z�HD�` AX=�B/�B+;�    ?z�HD�� AV��B.��B*�    ?z�HD�� AUr�B-q�B) *    ?z�HD�  ATiB,N�B'�Q    ?z�HD�` AR�JB+-B&�v    ?z�HD�� AQO�B*�B%��    ?z�HD�� AO�B(�.B$�v    ?z�HD�  AN�8B'�gB#��    ?z�HD�` AM?PB&��B"xC    ?z�HD�� AK�$B%��B!f3    ?z�HD�� AJ��B$��B X�    ?z�HD�  AI?�B#kuBJ#    ?z�HD�` AG�sB"XB<-    ?z�HD�� AF�.B!D�B1�    ?z�HD�� AEUQB 3B*     ?z�HD�  AD�B"�B%�    ?z�HD�` AB�zB�B#
    ?z�HD  AA��B�B|    ?z�HD�� A@K�B��B�    ?z�HD�  A?B�aB�    ?z�HD�` A=�tB�B#�    ?z�HDà A<��B��B&     ?z�HD�� A;Y�B�!B)�    ?z�HD�  A: RB�B1N    ?z�HD�` A8�CB��B6�    ?z�HDĠ A7��BڜB@i    ?z�HD�� A6{BۺBJ�    ?z�HD�  A5L�B�kBU�    ?z�HD�` A4�B�Bb3    ?z�HDŠ A2�B�Bos    ?z�HD�� A1�,B��B�    ?z�HD�  A0�}B��B�o    ?z�HD�` A/f�BxB
�?    ?z�HDƠ A.=�BB	�{    ?z�HD�� A-�B�B�!    ?z�HD�  A+��B+B�n    ?z�HD�` A*�]B
:IB��    ?z�HDǠ A)��B	LRB�    ?z�HD�� A(��B_IB)}    ?z�HD�  A'g(Bs+BD�    ?z�HD�` A&G�B�EB_�    ?z�HDȠ A%)�B�jB|�    ?z�HD�� A${B�XB�d    ?z�HD�  A"��BӣB ��    ?z�HD�` A!��B�A���    ?z�HDɠ A �xB'A���    ?z�HD�� A�IB*_A�F�    ?z�HD�  A�<B J,A���    ?z�HD�` A��A��A�ݸ    ?z�HDʠ At�A�eA�,4    ?z�HD�� Ae:A�b�A��    ?z�HD�  AWhA��
A�ӊ    ?z�HD�` AKKA���A�*;    ?z�HDˠ A@�A�J�A���    ?z�HD�� A4�A��A��X    ?z�HD�  A&�A��A�8s    ?z�HD�` A�A�HEA��    ?z�HD̠ AbA�/A��    ?z�HD�� AsA� A�Q    ?z�HD�  A CA�]�A�N    ?z�HD�` A��A��^A��    ?z�HD͠ A��A�%&Aㆅ    ?z�HD�� A�QA獥A��+    ?z�HD�  A�A��&A�^    ?z�HD�` A�{A�a�A��5    ?z�HDΠ A��A�ИA�A9    ?z�HD�� A
��A�@�A۳�    ?z�HD�  A
 DAߴyA�,�    ?z�HD�` A	5A�+#Aب<    ?z�HDϠ A�Aܤ�A�#�    ?z�HD�� A�A��Aբ�    ?z�HD�  A%YAٝIA�#�    ?z�HD�` A2�A�cAҧ�    ?z�HDР AAA֟A�.
    ?z�HD�� AR	A�#�A϶�    ?z�HD�  Ae�AӪHA�Dg    ?z�HD�` AzeA�3�A��X    ?z�HDѠ A ��Aо�A�b�    ?z�HD�� @�P�A�M3A���    ?z�HD�  @��fA�ޛAȌ�    ?z�HD�` @���A�pA�#'    ?z�HDҠ @���A�%Ażm    ?z�HD�� @�)KAɝ�A�Xx    ?z�HD�  @�e�A�7 A��c    ?z�HD�` @��.A���A��8    ?z�HDӠ @���A�p3A�8�    ?z�HD�� @�+�A��A���    ?z�HD�  @�r�A³�A��Y    ?z�HD�` @��-A�Y�A�,)    ?z�HDԠ @�\A�A��r    ?z�HD�� @�W�A��!A��+    ?z�HD�  @��A�VJA�43    ?z�HD�` @���A�A��    ?z�HDՠ @�R+A��RA��I    ?z�HD�� @��A�d2A�J�    ?z�HD�  @��A�7A� �    ?z�HD�` @�b�A��cA��    ?z�HD֠ @��#A���A�v    ?z�HD�� @�#�A�BrA�2�    ?z�HD�  @ۈ6A��VA��T    ?z�HD�` @��dA��_A���    ?z�HDנ @�VCA�}�A�sL    ?z�HD�� @��A�@*A�6�    ?z�HD�  @�-8A�~A���    ?z�HD�` @ӛ�A��eA���    ?z�HDؠ @��A��qA���    ?z�HD�� @�:A�^�A�WY    ?z�HD�  @��]A�-�A�%�    ?z�HD�` @�i�A��uA��    ?z�HD٠ @��cA���A��C    ?z�HD�� @�\A��KA���    ?z�HD�  @���A�s�A�g�    ?z�HD�` @�T�A�JkA�9�    ?z�HDڠ @��&A�"�A�/    ?z�HD�� @�T�A���A��    ?z�HD�  @��2A��VA��'    ?z�HD�` @�YKA���A���    ?z�HD۠ @��^A��}A�w�    ?z�HD�� @�bCA�{%A�S�    ?z�HD�  @��1A�^]A�0�    ?z�HD�` @�p�A�D�A��    ?z�HDܠ @��A�+A��    ?z�HD�� @��ZA�1A��a    ?z�HD�  @�1A���A��    ?z�HD�` @���A���A���    ?z�HDݠ @��A��/A�y�    ?z�HD�� @���A�͑A�]    ?z�HD�  @�2�A���A�A�    ?z�HD�` @��nA��{A�&&    ?z�HDޠ @�D�A���A�7    ?z�HD�� @��UA��CA���    ?z�HD�  @�Q�A���A���    ?z�HD�` @���A��UA��v    ?z�HDߠ @�XzA���A��&    ?z�HD�� @�ؽA���A�}�    ?z�HD�  @�T�A��>A�]�    ?z�HD�` @��pA���A�>�    ?z�HD� @�A[A���A�`    ?z�HD�� @��cA���A���    ?z�HD�  @�R�A��LA��    ?z�HD�` @�-?A��HA}�&    ?z�HD� @�
BA���A|,q    ?z�HD�� @���A���AzZ    ?z�HD�  @��BA�ͬAx�n    ?z�HD�` @��9A���Av�2    ?z�HD� @��A��_At�X    ?z�HD�� @�v�A��As4    ?z�HD�  @�^A~+�Aqq�    ?z�HD�` @�H�A|Y�Ao��    ?z�HD� @�4UAz�Am��    ?z�HD�� @�!�Ax��Al?_    ?z�HD�  @�-Av�-Aj��    ?z�HD�` @��Au)Ah�d    ?z�HD� @��Asc�Ag#/    ?z�HD�� @��=Aq��Aew?    ?z�HD�  @��oAo�Ac·    ?z�HD�` @���An'nAb&�    ?z�HD� @��XAln�A`��    ?z�HD�� @��Aj�A^�H    ?z�HD�  @���Ai@A]@X    ?z�HD�` @���AgX'A[�    ?z�HD� @��nAe�AZ
�    ?z�HD�� @�ђAd1AXs:    ?z�HD�  @���AbY#AV޸    ?z�HD�` @�ܓA`��AUM�    ?z�HD� @���A_�AS�    ?z�HD�� @���A]v4AR3�    ?z�HD�  @���A[��AP��    ?z�HD�` @��AZ@�AO$    ?z�HD� @��AX��AM��    ?z�HD�� @~IoAW�AL"�    ?z�HD�  @|l�AU�AJ��    ?z�HD�` @z�AS��AI&.    ?z�HD� @x��ARo�AG��    ?z�HD�� @v�AP��AF7A    ?z�HD�  @u
AOa�AD�    ?z�HD�` @sM�AMޛACP�    ?z�HD� @q��AL]#AA�I    ?z�HD�� @o��AJ�cA@s�    ?z�HD�  @m�^AIe5A?�    ?z�HD�` @l<�AG�,A=��    ?z�HD� @j~[AFu%A<<s    ?z�HD�� @h��AEgA:��    ?z�HD�  @gaAC�=A9|�    ?z�HD�` @eZ�AB#A87    ?z�HD� @c��A@��A6�    ?z�HD�� @a�DA?M�A5h    ?z�HD�  @`P�A=�	A4�    ?z�HD�` @^�wA<�uA2�"    ?z�HD�� @]�A;!tA1rw    ?z�HD�� @[��A9A01�    ?z�HD�  @Y��A8ePA.�S    ?z�HD�` @Xw�A7,A-�	    ?z�HD� @V�=A5��A,vK    ?z�HD�� @Us�A4]1A+<N    ?z�HD�  @S��A3
�A*P    ?z�HD�` @Rw�A1��A(��    ?z�HD� @P�DA0kA'�D    ?z�HD�� @O��A/A&d�    ?z�HD�  @N�A-�$A%4�    ?z�HD�` @L��A,��A$o    ?z�HD� @K(�A+FRA"�    ?z�HD�� @I�ZA*�A!�    ?z�HD�  @HK�A(´A �-    ?z�HD�` @F��A'��A]�    ?z�HD� @Ev�A&GA8U    ?z�HD�� @D�A%�A�    ?z�HD�  @B��A#��A�    ?z�HD�` @AI,A"��A��    ?z�HD� @?�A!k�A��    ?z�HD�� @>��A :�A��    ?z�HD�  @=.!A	�A}i    ?z�HD�` @;�AےAd�    ?z�HD� @:|A��AM    ?z�H
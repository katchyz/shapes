convert t_cell1_invitro_NAIN3.png t_cell2* t_oblong_* +append top.png
convert t_sphere_invitro_A.png c_* +append mid.png
convert t_sphere_invitro_C.png l_* +append bot.png
convert top.png mid.png bot.png -append start.png
rm top.png mid.png bot.png


convert t_cell1_invitro_NAIN3.png t_cell2*o.png t_oblong_*o.png +append top.png
convert t_sphere_invitro_A.png c_*o.png +append mid.png
convert t_sphere_invitro_C.png l_*o.png +append bot.png
convert top.png mid.png bot.png -append stop.png
rm top.png mid.png bot.png


convert t_cell1_invitro_NAIN3_TAG.png t_cell2*o_TAG.png t_oblong_*o_TAG.png +append top.png
convert t_sphere_invitro_A_TAG.png c_*o_TAG.png +append mid.png
convert t_sphere_invitro_C_TAG.png l_*o_TAG.png +append bot.png
convert top.png mid.png bot.png -append stop_TAG.png
rm top.png mid.png bot.png


convert tx_cell1_invitro_NAIN3_TREATED.png tx_cell2*_TREATED.png tx_oblong_*_TREATED.png +append top.png
convert tx_sphere_invitro_A_TREATED.png tx_*_CONTROL.png +append mid.png
convert tx_sphere_invitro_C_TREATED.png tx_*_LOG2RATIO.png +append bot.png
convert top.png mid.png bot.png -append utr5_cds_utr3.png
rm top.png mid.png bot.png

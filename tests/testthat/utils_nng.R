# ==============================================================================
# scclust for R -- R wrapper for the scclust library
# https://github.com/fsavje/scclust-R
#
# Copyright (C) 2016  Fredrik Savje -- http://fredriksavje.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/
# ==============================================================================

nng_clustering <- function(distance_object,
                           size_constraint,
                           seed_method = "exclusion_updating",
                           unassigned_method = "closest_seed",
                           radius = NULL,
                           primary_data_points = NULL,
                           secondary_unassigned_method = "ignore",
                           secondary_radius = NULL) {
  primary_radius <- "seed_radius"
  if (unassigned_method == "estimated_radius_closest_seed") {
    unassigned_method <- "closest_seed"
    primary_radius <- "estimated_radius"
  }

  if (secondary_unassigned_method == "estimated_radius_closest_seed") {
    secondary_unassigned_method <- "closest_seed"
    secondary_radius <- "estimated_radius"
  }

  make_clustering(distance_object,
                  size_constraint,
                  NULL,
                  NULL,
                  seed_method,
                  primary_data_points,
                  unassigned_method,
                  secondary_unassigned_method,
                  radius,
                  primary_radius,
                  secondary_radius,
                  NULL)
}


nng_clustering_batches <- function(distance_object,
                                   size_constraint,
                                   unassigned_method = "any_neighbor",
                                   radius = NULL,
                                   primary_data_points = NULL,
                                   batch_size = 100L) {
  primary_radius <- "seed_radius"
  if (unassigned_method == "estimated_radius_closest_seed") {
    unassigned_method <- "closest_seed"
    primary_radius <- "estimated_radius"
  }

  make_clustering(distance_object,
                  size_constraint,
                  NULL,
                  NULL,
                  "batches",
                  primary_data_points,
                  unassigned_method,
                  "ignore",
                  radius,
                  primary_radius,
                  NULL,
                  batch_size)
}


nng_clustering_types <- function(distance_object,
                                 type_labels,
                                 type_size_constraints,
                                 total_size_constraint = NULL,
                                 seed_method = "exclusion_updating",
                                 unassigned_method = "closest_seed",
                                 radius = NULL,
                                 primary_data_points = NULL,
                                 secondary_unassigned_method = "ignore",
                                 secondary_radius = NULL) {
  primary_radius <- "seed_radius"
  if (unassigned_method == "estimated_radius_closest_seed") {
    unassigned_method <- "closest_seed"
    primary_radius <- "estimated_radius"
  }

  if (secondary_unassigned_method == "estimated_radius_closest_seed") {
    secondary_unassigned_method <- "closest_seed"
    secondary_radius <- "estimated_radius"
  }

  make_clustering(distance_object,
                  total_size_constraint,
                  unclass(type_labels),
                  type_size_constraints,
                  seed_method,
                  primary_data_points,
                  unassigned_method,
                  secondary_unassigned_method,
                  radius,
                  primary_radius,
                  secondary_radius,
                  NULL)
}


test_nng_against_replica <- function(distance_object,
                                     size_constraint,
                                     seed_method,
                                     unassigned_method,
                                     radius,
                                     primary_data_points,
                                     secondary_unassigned_method,
                                     secondary_radius) {
  eval(bquote(expect_identical(nng_clustering(distance_object,
                                              size_constraint,
                                              seed_method,
                                              unassigned_method,
                                              radius,
                                              primary_data_points,
                                              secondary_unassigned_method,
                                              secondary_radius),
                               replica_nng_clustering(distance_object,
                                                      size_constraint,
                                                      seed_method,
                                                      unassigned_method,
                                                      radius,
                                                      primary_data_points,
                                                      secondary_unassigned_method,
                                                      secondary_radius))))
}


test_nng_batch_against_replica <- function(distance_object,
                                           size_constraint,
                                           unassigned_method,
                                           radius,
                                           primary_data_points,
                                           batch_size) {
  eval(bquote(expect_identical(nng_clustering_batches(distance_object,
                                                      size_constraint,
                                                      unassigned_method,
                                                      radius,
                                                      primary_data_points,
                                                      batch_size),
                               replica_nng_clustering_batches(distance_object,
                                                              size_constraint,
                                                              unassigned_method,
                                                              radius,
                                                              primary_data_points))))
}


test_nng_types_against_replica <- function(distance_object,
                                           type_labels,
                                           type_size_constraints,
                                           total_size_constraint,
                                           seed_method,
                                           unassigned_method,
                                           radius,
                                           primary_data_points,
                                           secondary_unassigned_method,
                                           secondary_radius) {
  eval(bquote(expect_identical(nng_clustering_types(distance_object,
                                                    type_labels,
                                                    type_size_constraints,
                                                    total_size_constraint,
                                                    seed_method,
                                                    unassigned_method,
                                                    radius,
                                                    primary_data_points,
                                                    secondary_unassigned_method,
                                                    secondary_radius),
                               replica_nng_clustering_types(distance_object,
                                                            type_labels,
                                                            type_size_constraints,
                                                            total_size_constraint,
                                                            seed_method,
                                                            unassigned_method,
                                                            radius,
                                                            primary_data_points,
                                                            secondary_unassigned_method,
                                                            secondary_radius))))
}


test_data <- matrix(c(9.8864788, 9.5334187, 1.4438035, 9.1983830, 0.2522823, 9.4571505, 6.5875638, 4.1881377, 8.8647318, 0.9832436,
                      2.5056677, 4.3667750, 7.7585977, 5.5608923, 1.4291665, 2.0975896, 4.7239460, 0.6701718, 0.5706785, 8.3357518,
                      8.8108044, 3.8732228, 1.2679993, 6.7329704, 2.7997385, 3.3545525, 7.5740223, 7.6224886, 5.5201533, 7.7733261,
                      9.4182215, 5.4427796, 8.5247363, 6.7010569, 7.8071688, 1.9711971, 8.3728024, 2.1925468, 2.5970255, 2.0737894,
                      0.2610044, 5.1401945, 0.4969963, 7.1052088, 8.2374700, 5.9078620, 7.0335312, 5.9402552, 5.3182458, 0.3577965,
                      4.1613895, 9.9921433, 4.7627655, 3.1090632, 1.8614876, 4.2361801, 6.4914996, 0.5896583, 1.1638259, 5.9283403,
                      8.7432867, 2.2935228, 9.4400448, 5.5108225, 1.4261551, 1.3681354, 8.0096757, 6.8969360, 7.2077437, 0.2096411,
                      6.3848557, 8.7883988, 3.6777667, 5.1399095, 9.5838920, 5.2293546, 9.3327463, 0.6138374, 8.2600290, 5.0601039,
                      7.9076553, 2.5775471, 0.1888883, 0.8374754, 9.7207443, 3.3814432, 2.1574054, 8.7780280, 8.8420260, 9.2293341,
                      2.6031073, 0.3215815, 7.9715905, 6.7020155, 0.5462608, 2.3495315, 4.6767730, 4.6767265, 5.3702017, 2.0973276,
                      4.7785926, 2.4512474, 4.5139801, 0.1272137, 1.6446873, 6.4288398, 4.2597485, 1.2558629, 8.3742285, 4.2803411,
                      4.0800889, 2.7986000, 2.4900506, 0.0709884, 6.3343779, 9.2610499, 5.7600991, 1.1279701, 3.6330739, 7.5678827,
                      8.3571291, 8.8736490, 4.5852803, 4.7344123, 6.8755419, 2.2833523, 9.1392932, 8.9908547, 0.2419371, 7.4014447,
                      3.0642209, 1.5105482, 0.5165190, 1.4992456, 9.2391617, 3.7435693, 9.6943672, 6.3859986, 1.4880940, 1.9614856,
                      2.1413810, 1.3862349, 8.4062852, 5.9949409, 0.4322822, 3.1052543, 6.8302958, 1.8398975, 2.5328603, 4.0329873,
                      2.9439773, 2.2481555, 7.0872746, 2.8510579, 2.1198794, 3.9406766, 0.1532272, 3.0352724, 5.2638885, 8.8737820,
                      1.8773351, 3.0337770, 1.0006110, 3.1454479, 8.7551967, 0.6299004, 7.7982707, 0.0601558, 1.3751565, 9.9118372,
                      2.1829523, 4.8971118, 2.2161192, 9.1081415, 9.7672505, 5.1161298, 5.8605465, 4.4583656, 4.2591806, 3.1056301,
                      5.6081774, 3.2256541, 9.4425947, 4.5824165, 4.6292410, 3.5696129, 3.1462250, 4.2502480, 4.8043871, 1.4110271,
                      8.5159693, 4.4613593, 7.8274567, 9.6075752, 1.5748424, 2.0887966, 4.3088750, 5.1206188, 8.5742208, 8.6563762,
                      1.1371549, 2.0677887, 0.8171307, 4.8362069, 5.8697862, 2.1282324, 5.2422277, 1.5011601, 0.2886923, 0.7483813,
                      8.7509141, 3.4195464, 6.0323489, 9.9495109, 3.1440682, 5.1336960, 3.7722057, 5.2393734, 7.7882531, 2.9442013,
                      8.1622170, 3.9695252, 5.2911191, 6.3593522, 8.7533614, 2.5513621, 0.0080694, 2.2518512, 7.1071369, 4.1024481,
                      5.2220798, 6.0162450, 6.8969638, 3.6450808, 0.9015877, 4.3675659, 2.1800387, 9.1784456, 7.9252138, 6.5599427,
                      4.1778268, 3.0190771, 4.6519740, 1.3599360, 0.7984144, 6.5314125, 6.7154715, 0.4307526, 7.0392524, 0.6159035,
                      2.7369246, 4.5515802, 0.3479291, 6.6874290, 4.1825261, 3.6355263, 1.4064719, 2.4076534, 4.5052295, 9.3771001,
                      9.9166656, 6.0561301, 2.5403798, 9.3045475, 2.7519763, 3.4406770, 7.4246680, 1.4024813, 8.4427075, 8.6061813,
                      7.1432499, 9.0203128, 8.2869790, 7.2126135, 9.3705428, 1.5285590, 4.3961239, 7.1097714, 4.1752013, 7.6998725,
                      6.6958014, 9.5274184, 2.5551244, 5.4565950, 1.0219692, 0.1885554, 4.4083655, 8.9724812, 1.9127612, 6.7924176,
                      9.0551638, 2.3861006, 2.6651538, 1.5010527, 1.0359159, 0.8056488, 5.4818822, 6.1313588, 2.4763447, 2.1361185,
                      3.2891006, 6.8067610, 1.9670535, 4.3007660, 9.4782948, 0.9102638, 4.1027709, 3.5208572, 4.1857452, 1.4308891,
                      0.0053061, 9.1491113, 9.2439385, 4.1658397, 4.5911336, 7.9839318, 4.9227881, 3.2086658, 9.8414663, 5.8867952,
                      9.2975727, 4.1343185, 8.7253190, 0.7000128, 4.6347384, 0.6665857, 2.3774150, 3.2917854, 5.7443710, 9.1875091,
                      0.8627268, 3.4166306, 4.6976863, 3.2177075, 7.2292536, 6.0964692, 3.9742966, 1.1350838, 8.6001734, 7.0424698,
                      0.2588156, 8.4971426, 3.3800580, 6.4639627, 7.7383366, 8.2269609, 7.3514310, 9.5345110, 8.5013446, 5.2166344,
                      9.8621898, 9.4074709, 0.7045638, 6.1124584, 5.6920918, 4.9846190, 1.7062121, 9.9636886, 2.0590992, 6.2415306,
                      4.8674707, 7.9304188, 7.6294340, 6.9386780, 4.0003315, 5.1619206, 7.1960504, 3.3643685, 0.2462535, 8.7758628,
                      7.3216265, 9.5352227, 2.8895670, 1.9376775, 8.2471345, 8.3593096, 0.7664363, 4.4944122, 2.7284463, 2.5643176,
                      5.9922859, 1.2448434, 6.9622285, 2.9651849, 6.8006760, 4.0043550, 3.0074378, 5.1175464, 3.5743238, 1.5580976,
                      3.6849451, 3.0377586, 1.0880876, 2.2381412, 3.2749316, 1.9172963, 3.0251030, 0.9169774, 2.3493001, 7.9258480,
                      8.9308971, 7.3642679, 5.2817366, 7.5958463, 8.3007121, 4.6393643, 6.2226338, 3.2002683, 6.7955531, 7.0025285,
                      2.5959244, 3.3793498, 8.8491551, 9.9342539, 3.9114550, 6.1787541, 1.2017288, 9.3212964, 3.8365175, 5.7830811,
                      4.5281600, 0.3005846, 3.0944057, 2.1364599, 6.9534450, 8.2692273, 4.2242151, 9.7132848, 8.2590850, 9.6823141,
                      7.8058451, 9.3412637, 8.3032893, 0.7173583, 4.3427734, 4.0839755, 7.0928779, 0.0965491, 9.7860159, 8.7683784,
                      5.7870079, 5.1424475, 1.0330792, 6.2974229, 6.4065907, 1.1486129, 1.2885714, 1.7048658, 6.2479521, 2.4399049,
                      2.8158086, 8.8743580, 7.7090446, 5.1876554, 0.4493401, 5.6795602, 1.0645745, 7.4346616, 7.2195078, 0.0870234,
                      7.2116999, 5.1253681, 2.6339692, 5.7229997, 3.2603872, 8.6472702, 5.0064336, 6.7279994, 3.9728728, 8.2376183,
                      1.2590379, 8.5980723, 2.4763922, 6.9293997, 9.3189098, 7.1329985, 1.2319888, 7.1014940, 6.6878973, 0.2905650,
                      7.6084892, 4.1940005, 8.9168762, 5.4901281, 5.3660451, 6.0819846, 0.6877188, 1.5963480, 6.0695931, 8.5143512,
                      6.2266475, 0.1466320, 2.1184612, 0.8735507, 9.1391024, 1.2220005, 6.3694136, 7.9719589, 0.7217801, 9.5892975,
                      8.1790164, 4.2688940, 7.0153176, 3.5538809, 9.8802535, 8.8403277, 8.7201174, 8.0994120, 8.2535592, 6.5442346,
                      6.8230444, 6.1298893, 9.5683539, 6.0581474, 5.2037278, 9.3806834, 7.2443792, 3.1653017, 6.2563758, 1.3229629,
                      6.9911344, 6.2283780, 2.4169716, 3.8388245, 1.5896208, 5.7465669, 7.7747985, 3.0426382, 4.5733009, 2.7593206,
                      1.7157791, 5.4243480, 6.7308652, 0.2804719, 7.7104778, 2.2577177, 1.8124760, 2.5855444, 3.6685277, 9.8576158,
                      8.0783601, 5.0065396, 3.3774424, 0.5406812, 0.3127767, 2.0566499, 1.5380621, 0.1784556, 5.5113529, 7.7177685,
                      6.0108807, 1.8774400, 5.2623750, 3.0403199, 3.2238705, 5.2395272, 1.8355189, 4.2771896, 2.7836595, 7.6104465,
                      5.5204927, 6.4776796, 6.7572047, 4.7796497, 5.7981008, 9.7522182, 3.5208639, 8.4369729, 4.4531533, 9.9057528,
                      7.1119885, 5.5543091, 3.8313569, 9.2930890, 6.0771053, 1.0263923, 4.9122370, 8.4121624, 2.8947766, 6.9664438,
                      0.7025950, 7.5574560, 3.4977492, 4.2776034, 3.8553654, 7.2755186, 4.2680234, 1.7741419, 6.3079911, 1.5032623,
                      2.6034652, 8.1788513, 8.3798712, 9.7665847, 6.7709674, 9.5074545, 4.0528027, 6.4577938, 6.9921386, 4.5514695,
                      7.9492246, 9.4643046, 2.1532076, 0.7690384, 1.7560918, 3.2743514, 5.4085739, 6.2442380, 0.8132703, 5.6105268,
                      7.7060267, 6.8598246, 3.7696891, 9.3217996, 9.3456232, 1.0153674, 0.4742949, 7.7425117, 4.9868881, 1.7707695,
                      4.9690850, 2.5150359, 1.7463894, 4.4097088, 4.7807428, 6.5742406, 4.7349653, 9.0545933, 0.3729086, 2.8960091,
                      0.6660246, 7.1273605, 5.2633905, 1.7367564, 4.0366838, 8.5842167, 8.2909883, 3.2793649, 6.3071343, 5.6778813,
                      9.5484407, 9.9923540, 8.2705353, 8.2225506, 4.6813513, 2.4829560, 5.3522673, 0.2050281, 5.7242142, 6.0597425,
                      3.0220492, 2.9161371, 1.5690025, 5.4516403, 2.1657126, 9.1295620, 9.4378111, 2.3629466, 4.0763082, 6.3448404,
                      9.3102445, 8.8249738, 0.8296446, 1.3150271, 3.2782584, 2.1409859, 2.0084680, 6.2864118, 3.4879060, 6.6926535,
                      9.1358863, 4.0561564, 6.2582506, 6.1049613, 3.8826148, 6.1550616, 9.3029332, 9.9275625, 4.4005931, 9.7800333,
                      5.6540450, 0.5673786, 5.4281337, 5.6133391, 2.9848171, 5.6206997, 4.1653964, 4.6621242, 9.5610451, 4.5465568,
                      1.5588946, 6.1168857, 0.8911647, 4.6820360, 3.8684992, 4.0652062, 5.1827969, 3.5580135, 4.5982369, 8.6401932,
                      6.0720160, 7.6055375, 3.7119019, 2.9681026, 9.4014829, 4.9748604, 0.1754220, 1.3548397, 1.9949136, 3.9147724,
                      6.2297976, 7.8373824, 3.6851039, 7.6699210, 8.3947480, 6.2201200, 3.0793077, 2.0629595, 0.9205899, 3.6312067,
                      1.6765476, 8.7917506, 8.4265479, 2.6722686, 2.8020531, 1.3892732, 4.3461764, 9.7622859, 0.2252866, 1.1286058,
                      8.5718858, 1.6273478, 1.1651023, 9.3776233, 4.7166524, 9.1281146, 9.2973019, 9.1233737, 2.8360941, 7.7251178,
                      7.9015563, 6.2096598, 7.4709683, 2.0535294, 2.8435980, 4.0612869, 5.3970550, 7.3284114, 5.2569261, 2.4069560,
                      9.0820023, 0.0179511, 6.9756183, 2.5234871, 5.1635663, 5.2386350, 0.5474039, 7.2367149, 9.7811719, 5.3821483,
                      2.5698939, 5.9813355, 4.5853144, 8.5008364, 4.6218968, 0.7272020, 4.6893724, 2.0475365, 3.0613638, 3.0905688,
                      9.4121383, 3.6418925, 7.2280507, 4.6160563, 2.1240530, 6.1078700, 1.6453180, 8.4615377, 4.5315998, 1.5850040,
                      6.1312230, 9.1588629, 6.8618556, 0.7439660, 4.6078704, 1.5750012, 6.8344732, 1.4818084, 1.0448251, 6.7359200,
                      0.6900269, 5.6731813, 9.9766540, 8.5046093, 6.1908063, 5.4683837, 2.4875887, 7.3739661, 5.5407294, 2.0234338,
                      5.8317572, 2.9294540, 3.3032952, 4.2809430, 9.5829027, 5.3207716, 9.8785553, 3.7571179, 8.5618081, 9.9607702,
                      4.1032397, 8.9964256, 2.1266050, 7.7473869, 2.2107254, 7.7086772, 8.1738754, 5.0176811, 8.8417945, 6.2222806,
                      8.3802380, 0.0466458, 7.0576358, 0.5311094, 8.1285884, 5.0759845, 3.2618509, 2.7266125, 0.0044663, 2.3860662,
                      7.2313710, 4.6025654, 8.0322270, 3.9407041, 2.3456653, 4.3426403, 8.4796601, 0.2388167, 7.0391442, 1.9171275,
                      4.9131876, 6.1380464, 7.8309760, 8.5391176, 6.4679985, 1.2487112, 3.0466185, 7.0228841, 6.8040858, 2.1037820,
                      7.5033218, 3.4477706, 4.7885792, 2.6627783, 5.3193467, 3.7458084, 1.3430758, 5.4886965, 0.4010457, 2.1731525,
                      2.0321975, 7.9641331, 0.3720254, 7.5784208, 5.0827124, 9.7990545, 3.8519725, 0.3346799, 4.4656122, 7.9364382,
                      2.9729036, 3.5321523, 3.4459069, 2.5596759, 8.6108970, 0.3868314, 9.2186846, 4.1528324, 2.7939028, 0.0069074,
                      6.7944193, 2.7184657, 8.4497151, 9.2311341, 4.5996791, 3.0903003, 7.9025319, 2.3677533, 8.5785200, 3.4631532,
                      0.4176288, 1.0881960, 0.3708095, 8.3898719, 7.8693205, 3.3418791, 1.8730164, 0.7059090, 3.6898606, 2.5659499,
                      6.0072536, 7.6580760, 7.0168410, 9.2831984, 4.2242990, 1.9963713, 9.3549982, 8.0175997, 1.3112327, 1.5688394,
                      4.2498301, 4.7339171, 2.5230025, 1.1977433, 8.3047549, 3.8421890, 0.0431153, 0.9043749, 1.2945730, 0.7282752,
                      7.2455059, 6.7604505, 0.8545563, 0.4253526, 5.5990128, 7.5269070, 5.2609363, 3.7796408, 9.9431848, 3.3668412,
                      6.1544839, 8.0625708, 3.7093271, 9.0770032, 4.9733893, 9.7593250, 8.1846358, 6.6496541, 0.9035978, 7.0473336,
                      6.7099209, 9.6892540, 3.0721679, 3.2680862, 2.4838190, 2.8453736, 1.4844134, 0.2325087, 4.9269801, 9.7255236,
                      3.1609354, 5.9504859, 6.4607226, 1.4588790, 4.0850679, 9.9494176, 3.7241753, 8.3946787, 6.3022761, 2.8481494,
                      1.4370508, 4.2719194, 4.8355073, 1.2761548, 8.7797514, 3.9266844, 5.8794236, 6.2284978, 1.4531219, 5.5755507,
                      5.6856727, 3.3051306, 0.5779532, 4.7102193, 5.6279902, 9.7107101, 7.3883444, 7.3343930, 6.6897779, 7.9945834,
                      4.7172250, 1.4076662, 3.4434990, 8.7564740, 9.8500421, 2.3985647, 9.6960007, 2.6405423, 1.9812322, 7.3725567,
                      2.8214478, 6.1441961, 1.3165908, 7.2738551, 4.9275615, 4.0848332, 1.4137695, 6.4601333, 8.1038091, 3.6685625), ncol = 2)

test_distances1 <- distances::distances(test_data)

primary_data_points <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
                         TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
                         FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE,
                         TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE,
                         TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
                         TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
                         FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
                         FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
                         TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE,
                         TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE,
                         FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
                         FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE,
                         FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE,
                         FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE,
                         TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE,
                         TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                         FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
                         FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
                         FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
                         TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
                         TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE,
                         TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE,
                         TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE,
                         TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE,
                         FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE,
                         FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
                         TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE,
                         FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
                         FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
                         TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
                         FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE,
                         FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
                         FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE,
                         FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE,
                         FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)

test_radius <- 0.2
type_test_radius <- 0.6

types1 <- as.integer(c(0, 1, 1, 0, 1, 2, 2, 0, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 1, 0, 0, 1, 0, 2, 3, 2, 1, 0, 2, 1, 1, 0, 2, 0,
                       1, 0, 1, 1, 1, 1, 2, 2, 1, 0, 1, 1, 1, 2, 1, 1, 1, 2, 0, 2, 1, 1, 2, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1,
                       1, 1, 2, 0, 0, 1, 1, 1, 2, 1, 1, 1, 3, 2, 2, 2, 0, 3, 3, 2, 2, 2, 1, 3, 2, 1, 3, 1, 2, 1, 2, 2, 1, 0,
                       0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 2, 1, 1, 2, 1, 2, 2, 1, 1, 1, 2, 1, 0, 2, 0, 1, 1, 3, 1, 1, 0, 0, 0,
                       2, 2, 1, 1, 0, 0, 1, 2, 1, 3, 0, 1, 1, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 0, 1, 1, 1, 1, 0, 2, 2, 2, 1, 0,
                       3, 0, 1, 0, 0, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 0, 2, 0, 2, 3, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2,
                       1, 2, 1, 2, 2, 2, 1, 0, 2, 1, 2, 2, 2, 0, 2, 0, 0, 3, 2, 2, 2, 1, 1, 2, 1, 1, 2, 3, 0, 0, 2, 0, 2, 2,
                       0, 1, 0, 0, 2, 0, 2, 1, 1, 1, 2, 0, 1, 1, 0, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 0, 2, 1, 2, 1, 2, 2,
                       0, 1, 1, 2, 1, 2, 2, 1, 3, 0, 2, 2, 0, 2, 1, 2, 2, 2, 1, 0, 3, 1, 1, 0, 1, 2, 1, 1, 0, 0, 2, 2, 1, 2,
                       1, 0, 2, 0, 1, 1, 0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 0, 0, 0, 1, 2, 1, 1, 2, 0, 1, 2, 0, 0, 1, 2, 3, 0, 1,
                       1, 1, 3, 0, 1, 0, 2, 1, 0, 1, 1, 2, 1, 2, 1, 0, 0, 2, 1, 2, 2, 1, 1, 1, 1, 1, 0, 1, 1, 2, 1, 0, 2, 2,
                       2, 2, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 2, 2, 2, 1, 0, 1, 2, 1, 1, 2, 1, 2, 2, 1,
                       1, 3, 1, 0, 2, 2, 0, 2, 1, 1, 1, 0, 2, 2, 1, 2, 1, 2, 1, 0, 2, 1, 1, 2, 1, 0, 2, 1, 0, 0, 1, 3, 0, 2,
                       0, 1, 3, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 2, 1, 0, 1, 3, 2, 1, 1, 1, 2, 0, 2, 1, 1, 1, 2, 2, 1, 2, 1,
                       1, 1, 1, 2, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 1, 3, 1, 1, 1, 1, 0, 1, 2, 1))

types2 <- as.integer(c(0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0,
                       0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0,
                       1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0,
                       0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1,
                       1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0,
                       1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0,
                       1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1,
                       0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0,
                       1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1,
                       1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0,
                       0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0,
                       1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
                       1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1,
                       1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0))

library(Rscclust)
context("greedy_clustering.R")

source("../replica/replica_greedy.R", local = TRUE)

sound_distance_obj <- make_distances(matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1), ncol = 1))
unsound_distance_obj <- structure(matrix(letters[1:9], ncol = 3),
                                  ids = NULL,
                                  normalization = diag(3),
                                  weights = diag(3),
                                  class = c("Rscc_distances"))
sound_size_constraint <- 3
unsound_size_constraint <- 0
sound_batch_assign <- TRUE
unsound_batch_assign <- "A"
sound_clustering <- Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a"))
unsound_clustering <- make_Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a"), 3, NULL)
sound_deep_copy <- TRUE
unsound_deep_copy <- "A"

test_that("greedy clustering functions check input.", {
  expect_match(class(top_down_greedy_clustering(sound_distance_obj,
                                                sound_size_constraint,
                                                sound_batch_assign,
                                                sound_clustering)), "Rscc_clustering", fixed = TRUE)
  expect_error(top_down_greedy_clustering(unsound_distance_obj,
                                          sound_size_constraint,
                                          sound_batch_assign,
                                          sound_clustering))
  expect_error(top_down_greedy_clustering(sound_distance_obj,
                                          unsound_size_constraint,
                                          sound_batch_assign,
                                          sound_clustering))
  expect_error(top_down_greedy_clustering(sound_distance_obj,
                                          sound_size_constraint,
                                          unsound_batch_assign,
                                          sound_clustering))
  expect_error(top_down_greedy_clustering(sound_distance_obj,
                                          sound_size_constraint,
                                          sound_batch_assign,
                                          unsound_clustering))

  expect_match(class(top_down_greedy_clustering_internal(sound_distance_obj,
                                                         sound_size_constraint,
                                                         sound_batch_assign,
                                                         sound_clustering,
                                                         sound_deep_copy)), "Rscc_clustering", fixed = TRUE)
  expect_error(top_down_greedy_clustering_internal(unsound_distance_obj,
                                                   sound_size_constraint,
                                                   sound_batch_assign,
                                                   sound_clustering,
                                                   sound_deep_copy))
  expect_error(top_down_greedy_clustering_internal(sound_distance_obj,
                                                   unsound_size_constraint,
                                                   sound_batch_assign,
                                                   sound_clustering,
                                                   sound_deep_copy))
  expect_error(top_down_greedy_clustering_internal(sound_distance_obj,
                                                   sound_size_constraint,
                                                   unsound_batch_assign,
                                                   sound_clustering,
                                                   sound_deep_copy))
  expect_error(top_down_greedy_clustering_internal(sound_distance_obj,
                                                   sound_size_constraint,
                                                   sound_batch_assign,
                                                   unsound_clustering,
                                                   sound_deep_copy))
  expect_error(top_down_greedy_clustering_internal(sound_distance_obj,
                                                   sound_size_constraint,
                                                   sound_batch_assign,
                                                   sound_clustering,
                                                   unsound_deep_copy))
})

test_data <- matrix(c(0.0436, 0.9723, 0.5366, 0.1065, 0.5340, 0.3437, 0.2933, 0.4599, 0.2895, 0.2217,
                      0.9043, 0.9513, 0.6091, 0.5963, 0.0520, 0.1248, 0.7416, 0.4801, 0.4345, 0.5842,
                      0.1488, 0.2255, 0.1243, 0.5040, 0.6102, 0.4197, 0.4248, 0.1986, 0.2409, 0.7515,
                      0.3019, 0.2094, 0.5941, 0.7086, 0.3111, 0.7290, 0.3675, 0.3241, 0.5645, 0.1056,
                      0.5747, 0.5449, 0.3050, 0.1893, 0.9559, 0.4274, 0.1691, 0.3504, 0.5032, 0.4387,
                      0.1264, 0.5995, 0.4948, 0.8134, 0.1887, 0.4399, 0.2289, 0.4345, 0.7754, 0.5619,
                      0.3135, 0.0413, 0.5932, 0.4756, 0.3365, 0.7892, 0.0705, 0.2215, 0.8149, 0.6639,
                      0.2179, 0.6081, 0.1530, 0.7410, 0.1611, 0.1994, 0.1910, 0.9720, 0.3889, 0.9728,
                      0.9022, 0.2899, 0.7509, 0.5702, 0.1775, 0.9618, 0.5706, 0.1021, 0.6040, 0.5551,
                      0.2679, 0.1904, 0.9281, 0.8282, 0.2820, 0.5084, 0.9033, 0.8978, 0.4608, 0.8428,
                      0.7790, 0.1411, 0.2035, 0.8918, 0.7072, 0.5913, 0.3998, 0.1688, 0.1406, 0.6926,
                      0.1030, 0.2569, 0.9865, 0.5749, 0.5187, 0.7032, 0.1908, 0.3070, 0.2905, 0.0458,
                      0.4331, 0.5083, 0.7734, 0.3279, 0.4797, 0.2472, 0.5653, 0.4473, 0.2393, 0.2677,
                      0.7235, 0.9135, 0.8289, 0.5878, 0.0628, 0.5516, 0.0953, 0.7469, 0.2721, 0.5620,
                      0.7917, 0.0066, 0.4844, 0.1268, 0.9351, 0.5901, 0.3497, 0.5899, 0.1751, 0.1245,
                      0.2695, 0.2760, 0.2081, 0.2762, 0.1797, 0.0863, 0.8187, 0.9765, 0.4429, 0.7980,
                      0.5352, 0.2207, 0.9821, 0.3811, 0.0052, 0.9403, 0.3465, 0.8163, 0.3416, 0.1401,
                      0.7982, 0.5720, 0.0764, 0.3707, 0.2718, 0.2511, 0.7885, 0.9428, 0.7657, 0.9053,
                      0.2173, 0.2734, 0.0150, 0.2673, 0.0675, 0.7435, 0.2608, 0.8335, 0.4047, 0.0449,
                      0.3815, 0.2003, 0.3374, 0.3810, 0.3093, 0.8579, 0.4245, 0.5590, 0.7313, 0.1392,
                      0.7651, 0.0990, 0.8497, 0.6004, 0.4276, 0.0759, 0.8396, 0.6570, 0.8217, 0.3836,
                      0.2844, 0.4514, 0.8184, 0.8845, 0.9190, 0.6074, 0.5558, 0.2583, 0.0977, 0.1340,
                      0.7432, 0.6822, 0.2847, 0.4911, 0.8880, 0.1184, 0.4989, 0.8032, 0.6957, 0.6037,
                      0.9236, 0.0298, 0.9962, 0.2664, 0.3726, 0.3677, 0.7628, 0.8168, 0.0561, 0.4812,
                      0.7474, 0.0942, 0.4603, 0.5352, 0.1658, 0.9764, 0.2241, 0.5676, 0.3937, 0.0339,
                      0.6563, 0.0607, 0.4649, 0.2805, 0.3180, 0.9205, 0.4763, 0.5221, 0.8002, 0.7053,
                      0.0569, 0.3173, 0.6423, 0.6311, 0.5560, 0.2583, 0.9154, 0.6428, 0.4902, 0.9706,
                      0.8906, 0.4008, 0.0853, 0.9762, 0.7717, 0.2763, 0.2713, 0.2722, 0.0866, 0.0508,
                      0.3016, 0.4123, 0.2925, 0.6964, 0.1037, 0.0330, 0.8935, 0.8227, 0.5834, 0.7269,
                      0.6629, 0.8114, 0.7532, 0.5804, 0.7616, 0.2694, 0.0428, 0.7617, 0.8282, 0.1830,
                      0.2099, 0.6619, 0.8957, 0.6263, 0.4265, 0.4317, 0.1860, 0.2142, 0.6227, 0.8592,
                      0.6843, 0.3839, 0.7468, 0.0516, 0.6644, 0.5902, 0.9871, 0.3640, 0.9319, 0.3305,
                      0.6103, 0.8239, 0.2603, 0.7243, 0.7425, 0.1016, 0.2734, 0.7640, 0.6924, 0.7348,
                      0.8800, 0.0876, 0.0821, 0.8307, 0.3544, 0.3124, 0.7081, 0.9131, 0.0557, 0.2447,
                      0.3204, 0.1114, 0.1015, 0.9235, 0.3113, 0.5089, 0.7095, 0.7806, 0.8769, 0.6782,
                      0.7728, 0.9863, 0.5632, 0.3743, 0.7134, 0.7039, 0.6481, 0.6429, 0.7933, 0.6501,
                      0.9110, 0.1240, 0.1708, 0.7673, 0.4043, 0.8887, 0.4868, 0.2946, 0.9142, 0.5394,
                      0.9479, 0.3096, 0.2727, 0.9137, 0.8863, 0.9460, 0.9094, 0.5932, 0.1355, 0.7493,
                      0.5467, 0.6732, 0.0740, 0.7452, 0.0681, 0.4331, 0.4355, 0.7166, 0.0614, 0.4220,
                      0.0303, 0.5097, 0.7536, 0.6582, 0.1357, 0.9994, 0.2254, 0.9082, 0.9788, 0.6023,
                      0.4076, 0.9632, 0.7511, 0.5429, 0.1457, 0.5213, 0.5214, 0.9580, 0.2671, 0.7701,
                      0.0094, 0.4813, 0.3005, 0.1531, 0.1656, 0.9927, 0.5729, 0.6650, 0.5619, 0.9005,
                      0.8632, 0.9320, 0.8907, 0.9360, 0.7495, 0.7828, 0.1481, 0.3862, 0.0058, 0.3169,
                      0.1705, 0.5135, 0.6803, 0.0532, 0.3539, 0.0576, 0.1468, 0.6572, 0.4748, 0.9416,
                      0.6154, 0.0761, 0.5800, 0.8364, 0.8972, 0.6302, 0.9241, 0.7039, 0.0604, 0.9804,
                      0.6092, 0.6828, 0.2491, 0.3101, 0.4305, 0.1709, 0.8160, 0.2603, 0.2984, 0.1309,
                      0.3241, 0.1859, 0.4224, 0.1646, 0.6320, 0.4228, 0.5165, 0.4340, 0.3184, 0.0033,
                      0.4847, 0.4590, 0.3563, 0.2993, 0.8750, 0.1273, 0.8160, 0.4488, 0.2434, 0.8626,
                      0.8613, 0.7686, 0.8805, 0.4834, 0.0311, 0.4737, 0.5532, 0.1091, 0.7885, 0.1027,
                      0.7573, 0.4336, 0.5853, 0.9770, 0.3554, 0.1929, 0.7408, 0.5537, 0.3899, 0.6372,
                      0.9506, 0.7599, 0.9998, 0.6086, 0.7185, 0.4918, 0.2772, 0.3168, 0.2936, 0.2721,
                      0.3825, 0.2715, 0.6898, 0.6994, 0.7635, 0.3311, 0.5181, 0.4884, 0.3868, 0.6931,
                      0.1077, 0.3666, 0.3665, 0.1460, 0.8484, 0.4961, 0.1908, 0.4340, 0.5253, 0.6901,
                      0.7609, 0.8985, 0.6216, 0.2682, 0.6609, 0.2153, 0.9731, 0.0020, 0.3362, 0.6747,
                      0.3465, 0.8292, 0.9293, 0.5028, 0.6093, 0.7261, 0.3670, 0.4184, 0.8262, 0.2517,
                      0.0499, 0.3679, 0.4195, 0.9439, 0.1771, 0.0210, 0.3202, 0.7672, 0.2591, 0.9103,
                      0.8023, 0.5417, 0.3385, 0.7690, 0.7808, 0.4472, 0.5012, 0.9965, 0.3661, 0.3693,
                      0.3317, 0.2329, 0.9374, 0.3738, 0.3931, 0.9351, 0.5914, 0.0959, 0.6950, 0.6815,
                      0.2915, 0.3797, 0.2768, 0.9789, 0.3148, 0.3074, 0.4861, 0.5774, 0.4765, 0.2996,
                      0.3532, 0.7821, 0.5237, 0.9317, 0.6447, 0.5887, 0.9068, 0.0493, 0.6286, 0.1366,
                      0.7414, 0.5123, 0.3618, 0.1499, 0.1227, 0.3585, 0.5319, 0.1420, 0.4395, 0.2862,
                      0.6642, 0.8716, 0.4487, 0.1079, 0.0287, 0.7401, 0.3932, 0.8855, 0.0883, 0.2991,
                      0.5674, 0.8675, 0.7577, 0.7437, 0.0866, 0.1912, 0.8964, 0.8285, 0.7130, 0.8011,
                      0.5124, 0.7860, 0.7177, 0.0468, 0.1683, 0.5912, 0.2501, 0.0570, 0.9274, 0.5138,
                      0.5329, 0.5956, 0.8746, 0.2677, 0.5931, 0.9333, 0.6614, 0.2535, 0.6653, 0.1488,
                      0.6792, 0.0936, 0.9801, 0.0703, 0.5279, 0.4895, 0.0050, 0.8447, 0.0952, 0.5688,
                      0.6538, 0.8613, 0.2261, 0.6304, 0.7880, 0.1550, 0.3924, 0.3985, 0.8525, 0.2062,
                      0.8634, 0.2468, 0.0525, 0.5784, 0.8691, 0.7208, 0.7920, 0.2212, 0.1496, 0.4202,
                      0.0469, 0.7316, 0.9837, 0.5906, 0.0382, 0.5070, 0.8849, 0.6483, 0.2662, 0.4067,
                      0.0086, 0.2154, 0.8864, 0.2034, 0.1259, 0.0728, 0.6471, 0.3750, 0.1946, 0.0092,
                      0.4901, 0.1939, 0.0974, 0.3322, 0.2288, 0.2620, 0.2512, 0.2771, 0.9667, 0.9129,
                      0.2792, 0.7201, 0.4196, 0.6935, 0.0872, 0.9529, 0.2404, 0.4489, 0.6862, 0.3641,
                      0.4323, 0.1591, 0.5840, 0.5834, 0.6025, 0.1555, 0.8824, 0.9232, 0.5303, 0.8451,
                      0.5560, 0.2730, 0.1702, 0.6730, 0.1005, 0.2427, 0.5696, 0.1472, 0.9614, 0.6551,
                      0.4816, 0.8049, 0.5184, 0.6602, 0.8632, 0.4133, 0.7895, 0.7826, 0.2318, 0.2682,
                      0.0571, 0.0230, 0.2558, 0.8047, 0.2776, 0.0103, 0.9714, 0.6285, 0.0202, 0.8571,
                      0.6637, 0.6044, 0.6710, 0.4321, 0.4363, 0.3492, 0.8818, 0.7203, 0.8245, 0.5944,
                      0.3636, 0.0143, 0.4956, 0.6719, 0.5838, 0.8819, 0.5894, 0.6857, 0.9830, 0.1190,
                      0.9809, 0.1861, 0.4821, 0.1085, 0.7325, 0.5013, 0.9576, 0.0318, 0.2253, 0.4908,
                      0.2108, 0.3915, 0.7860, 0.3720, 0.5158, 0.4369, 0.3297, 0.1730, 0.7970, 0.0723,
                      0.8408, 0.3938, 0.4611, 0.3007, 0.1794, 0.5379, 0.7658, 0.5571, 0.8045, 0.0475,
                      0.6559, 0.3278, 0.3920, 0.7821, 0.4875, 0.6378, 0.9600, 0.1944, 0.8123, 0.5501,
                      0.7793, 0.8661, 0.3155, 0.4735, 0.8022, 0.8871, 0.7778, 0.6761, 0.6833, 0.4577,
                      0.5250, 0.0422, 0.5466, 0.2655, 0.9254, 0.3755, 0.5933, 0.9638, 0.1785, 0.5856,
                      0.9349, 0.0683, 0.4285, 0.5025, 0.5703, 0.2405, 0.3346, 0.2267, 0.4349, 0.7568,
                      0.5310, 0.9895, 0.0372, 0.3669, 0.1851, 0.4993, 0.8388, 0.3827, 0.9331, 0.1162,
                      0.6130, 0.2976, 0.9103, 0.2054, 0.8658, 0.6464, 0.1606, 0.1596, 0.3098, 0.8589,
                      0.1489, 0.2684, 0.7811, 0.3992, 0.0461, 0.1663, 0.2696, 0.0253, 0.8906, 0.0161,
                      0.9709, 0.2527, 0.2579, 0.2572, 0.5221, 0.8990, 0.0138, 0.5799, 0.8787, 0.5140,
                      0.6201, 0.4744, 0.1558, 0.5651, 0.4509, 0.4345, 0.7643, 0.7983, 0.7472, 0.7787,
                      0.0219, 0.3235, 0.3356, 0.7447, 0.2495, 0.6689, 0.8928, 0.9761, 0.1129, 0.1161,
                      0.2153, 0.6733, 0.4567, 0.1196, 0.0717, 0.5004, 0.2219, 0.9613, 0.1528, 0.1672,
                      0.1543, 0.4521, 0.9111, 0.2158, 0.5793, 0.8378, 0.6415, 0.8200, 0.9119, 0.0964,
                      0.7444, 0.7829, 0.0130, 0.7598, 0.1832, 0.9318, 0.5191, 0.5890, 0.3526, 0.6591,
                      0.9645, 0.0609, 0.9011, 0.9993, 0.9311, 0.4316, 0.4345, 0.9373, 0.5544, 0.9956,
                      0.9623, 0.0552, 0.1168, 0.5721, 0.6144, 0.3810, 0.5661, 0.8631, 0.0519, 0.3982,
                      0.0606, 0.4301, 0.9740, 0.1933, 0.8117, 0.3697, 0.3979, 0.9366, 0.6409, 0.7621,
                      0.1431, 0.8020, 0.0426, 0.2259, 0.3280, 0.5776, 0.6302, 0.4370, 0.6054, 0.1274,
                      0.6781, 0.6921, 0.4024, 0.0273, 0.3808, 0.2799, 0.2430, 0.1689, 0.1598, 0.0208,
                      0.5561, 0.8567, 0.0660, 0.4726, 0.5772, 0.0277, 0.0141, 0.1749, 0.4614, 0.9239), ncol = 2)

test_distances1 <- make_distances(test_data)
prev_clust1 <- Rscc_clustering(c(9, 3, 4, 4, 7, 0, 5, 0, 6, 7, 4, 7, 4, 6, 4, 9, 9, 2, 4, 6, 7, 5, 5, 5, 3, 5, 2, 9, 6, 4, 6, 8, 0,
                                 9, 2, 2, 5, 5, 5, 6, 5, 2, 1, 1, 3, 6, 7, 9, 1, 8, 4, 8, 3, 1, 6, 8, 5, 7, 0, 7, 1, 1, 3, 0, 3, 9,
                                 0, 9, 9, 2, 4, 9, 8, 2, 0, 4, 9, 7, 9, 9, 2, 1, 7, 9, 8, 4, 9, 4, 6, 8, 2, 6, 8, 8, 2, 8, 9, 9, 5,
                                 9, 4, 4, 8, 8, 6, 5, 5, 7, 8, 7, 2, 5, 4, 5, 2, 2, 4, 6, 7, 8, 2, 7, 5, 5, 6, 7, 2, 4, 9, 1, 6, 8,
                                 8, 1, 4, 9, 2, 7, 6, 1, 1, 5, 5, 7, 5, 2, 9, 8, 2, 1, 6, 7, 0, 6, 6, 8, 8, 0, 3, 7, 6, 0, 1, 4, 3,
                                 3, 1, 3, 8, 3, 0, 1, 6, 5, 9, 1, 0, 3, 9, 1, 6, 2, 3, 9, 5, 2, 9, 3, 6, 6, 2, 4, 8, 3, 5, 4, 1, 0,
                                 3, 2, 1, 4, 6, 4, 7, 3, 3, 7, 0, 8, 7, 3, 1, 2, 4, 3, 7, 4, 0, 6, 7, 9, 0, 0, 0, 7, 5, 8, 1, 9, 6,
                                 6, 6, 1, 3, 2, 0, 5, 1, 9, 1, 9, 1, 1, 3, 7, 1, 9, 1, 9, 8, 1, 2, 8, 8, 8, 0, 5, 5, 2, 6, 3, 2, 3,
                                 8, 2, 3, 2, 0, 7, 9, 7, 1, 0, 5, 1, 1, 1, 1, 7, 4, 2, 4, 9, 9, 6, 5, 8, 9, 2, 8, 2, 8, 3, 6, 7, 3,
                                 8, 6, 1, 7, 7, 6, 1, 5, 8, 6, 6, 6, 3, 7, 6, 5, 3, 7, 6, 7, 4, 4, 2, 1, 3, 0, 5, 3, 7, 5, 1, 0, 0,
                                 3, 7, 0, 1, 3, 7, 6, 6, 8, 9, 1, 5, 3, 6, 9, 9, 2, 7, 2, 9, 6, 8, 4, 2, 0, 1, 8, 1, 1, 3, 1, 7, 2,
                                 3, 4, 9, 8, 3, 6, 3, 9, 1, 6, 4, 3, 1, 8, 0, 2, 8, 6, 5, 2, 9, 0, 8, 7, 9, 9, 5, 7, 4, 2, 4, 0, 6,
                                 5, 3, 9, 5, 4, 1, 8, 4, 2, 6, 1, 1, 3, 8, 1, 6, 7, 9, 1, 0, 4, 3, 3, 7, 9, 4, 2, 0, 2, 7, 6, 0, 7,
                                 3, 0, 5, 9, 9, 6, 6, 8, 3, 1, 5, 3, 5, 6, 5, 1, 3, 6, 5, 2, 4, 5, 3, 2, 1, 6, 5, 3, 0, 3, 3, 8, 0,
                                 4, 0, 9, 3, 6, 1, 1, 1, 9, 5, 0, 3, 8, 5, 4, 6, 7, 8, 6, 5, 1, 4, 7, 3, 0, 0, 0, 7, 4, 6, 5, 1, 3,
                                 8, 0, 5, 4, 8))
prev_clust2 <- Rscc_clustering(c(17, 6, 5, 6, 19, 2, 8, 12, 3, 3, 14, 16, 1, 1, 7, 7, 4, 5, 0, 13, 19, 5, 9, 17, 8, 6, 8, 3, 8, 18, 15,
                                 14, 11, 5, 12, 17, 8, 4, 11, 13, 8, 0, 17, 17, 17, 5, 12, 15, 2, 13, 3, 13, 1, 5, 12, 9, 19, 0, 19, 16,
                                 19, 5, 2, 4, 17, 13, 9, 15, 1, 17, 13, 14, 11, 8, 11, 15, 14, 19, 15, 15, 4, 13, 2, 15, 3, 6, 4, 12, 7,
                                 14, 19, 6, 9, 11, 7, 7, 11, 10, 15, 10, 17, 5, 6, 0, 11, 10, 15, 4, 8, 15, 12, 4, 13, 0, 2, 17, 17, 0,
                                 13, 13, 7, 9, 3, 8, 19, 1, 3, 18, 12, 6, 5, 3, 14, 0, 19, 2, 0, 17, 12, 12, 17, 14, 2, 11, 17, 7, 14,
                                 1, 1, 12, 16, 3, 6, 19, 6, 2, 8, 17, 2, 17, 14, 14, 18, 9, 17, 12, 4, 14, 16, 12, 12, 19, 13, 7, 10, 18,
                                 16, 14, 11, 5, 2, 10, 15, 17, 3, 15, 16, 16, 3, 11, 6, 6, 16, 0, 13, 7, 5, 15, 14, 8, 11, 9, 14, 7, 7,
                                 15, 13, 6, 10, 4, 14, 16, 16, 7, 18, 0, 9, 11, 7, 18, 9, 17, 1, 5, 4, 12, 7, 1, 5, 10, 13, 17, 17, 13,
                                 8, 15, 0, 13, 0, 6, 4, 18, 14, 1, 4, 10, 1, 15, 2, 9, 1, 16, 7, 3, 3, 11, 7, 3, 4, 15, 9, 6, 17, 14, 4,
                                 16, 8, 14, 11, 2, 12, 16, 6, 9, 12, 19, 19, 1, 4, 15, 7, 2, 6, 0, 19, 19, 13, 16, 2, 4, 11, 10, 16, 8,
                                 19, 12, 15, 3, 11, 12, 11, 10, 0, 15, 16, 13, 12, 4, 3, 6, 16, 19, 3, 18, 9, 1, 13, 17, 12, 7, 8, 14, 2,
                                 3, 7, 9, 7, 10, 3, 14, 9, 10, 13, 7, 6, 18, 12, 11, 2, 10, 12, 3, 9, 15, 5, 10, 3, 1, 8, 10, 10, 6, 6, 6,
                                 6, 15, 13, 2, 12, 0, 12, 13, 11, 6, 17, 1, 0, 16, 17, 10, 0, 10, 11, 7, 4, 0, 15, 0, 16, 5, 7, 2, 13, 11,
                                 6, 19, 5, 11, 3, 17, 18, 15, 6, 8, 9, 6, 0, 13, 15, 17, 1, 5, 11, 15, 6, 6, 2, 11, 19, 5, 3, 7, 0, 0, 1, 5,
                                 11, 12, 16, 12, 14, 9, 14, 16, 19, 17, 15, 8, 12, 5, 2, 1, 16, 7, 17, 6, 15, 8, 4, 3, 11, 3, 13, 1, 6, 18,
                                 4, 18, 12, 9, 0, 19, 10, 9, 18, 3, 13, 8, 7, 3, 16, 7, 0, 11, 15, 16, 17, 16, 16, 17, 11, 13, 14, 13, 7, 10,
                                 4, 0, 2, 18, 17, 10, 4, 12, 14, 9, 18, 16, 8, 16, 11, 2, 19, 12, 8, 11, 18, 3, 7, 2))

test_that("greedy clustering functions cluster correctly.", {
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 2L,
                                              batch_assign = TRUE),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 2L,
                                                      batch_assign = TRUE))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 3L,
                                              batch_assign = TRUE),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 3L,
                                                      batch_assign = TRUE))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 10L,
                                              batch_assign = TRUE),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 10L,
                                                      batch_assign = TRUE))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 11L,
                                              batch_assign = TRUE),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 11L,
                                                      batch_assign = TRUE))


  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 2L,
                                              batch_assign = FALSE),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 2L,
                                                      batch_assign = FALSE))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 3L,
                                              batch_assign = FALSE),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 3L,
                                                      batch_assign = FALSE))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 10L,
                                              batch_assign = FALSE),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 10L,
                                                      batch_assign = FALSE))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 11L,
                                              batch_assign = FALSE),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 11L,
                                                      batch_assign = FALSE))


  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 2L,
                                              batch_assign = TRUE,
                                              existing_clustering = prev_clust1),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 2L,
                                                      batch_assign = TRUE,
                                                      existing_clustering = prev_clust1))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 3L,
                                              batch_assign = TRUE,
                                              existing_clustering = prev_clust1),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 3L,
                                                      batch_assign = TRUE,
                                                      existing_clustering = prev_clust1))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 2L,
                                              batch_assign = FALSE,
                                              existing_clustering = prev_clust1),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 2L,
                                                      batch_assign = FALSE,
                                                      existing_clustering = prev_clust1))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 3L,
                                              batch_assign = FALSE,
                                              existing_clustering = prev_clust1),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 3L,
                                                      batch_assign = FALSE,
                                                      existing_clustering = prev_clust1))


  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 2L,
                                              batch_assign = TRUE,
                                              existing_clustering = prev_clust2),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 2L,
                                                      batch_assign = TRUE,
                                                      existing_clustering = prev_clust2))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 3L,
                                              batch_assign = TRUE,
                                              existing_clustering = prev_clust2),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 3L,
                                                      batch_assign = TRUE,
                                                      existing_clustering = prev_clust2))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 2L,
                                              batch_assign = FALSE,
                                              existing_clustering = prev_clust2),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 2L,
                                                      batch_assign = FALSE,
                                                      existing_clustering = prev_clust2))
  expect_identical(top_down_greedy_clustering(distance_object = test_distances1,
                                              size_constraint = 3L,
                                              batch_assign = FALSE,
                                              existing_clustering = prev_clust2),
                   replica_top_down_greedy_clustering(distance_object = test_distances1,
                                                      size_constraint = 3L,
                                                      batch_assign = FALSE,
                                                      existing_clustering = prev_clust2))
})


test_data <- matrix(c(0.0436, 0.9723, 0.5366, 0.1065, 0.5340, 0.3437, 0.2933, 0.4599, 0.2895, 0.2217,
                      0.9043, 0.9513, 0.6091, 0.5963, 0.0520, 0.1248, 0.7416, 0.4801, 0.4345, 0.5842,
                      0.1488, 0.2255, 0.1243, 0.5040, 0.6102, 0.4197, 0.4248, 0.1986, 0.2409, 0.7515,
                      0.3019, 0.2094, 0.5941, 0.7086, 0.3111, 0.7290, 0.3675, 0.3241, 0.5645, 0.1056,
                      0.5747, 0.5449, 0.3050, 0.1893, 0.9559, 0.4274, 0.1691, 0.3504, 0.5032, 0.4387,
                      0.1264, 0.5995, 0.4948, 0.8134, 0.1887, 0.4399, 0.2289, 0.4345, 0.7754, 0.5619), ncol = 2)

test_distances2 <- make_distances(test_data)
prev_clust3ref <- Rscc_clustering(c(2, 0, 3, 0, 2, 0, 2, 0, 2, 2, 1, 2, 3, 1, 0, 3, 0, 1, 1, 2, 3, 1, 1, 3, 2, 1, 2, 3, 0, 0))
prev_clust3a <- Rscc_clustering(c(2, 0, 3, 0, 2, 0, 2, 0, 2, 2, 1, 2, 3, 1, 0, 3, 0, 1, 1, 2, 3, 1, 1, 3, 2, 1, 2, 3, 0, 0))
prev_clust3b <- Rscc_clustering(c(2, 0, 3, 0, 2, 0, 2, 0, 2, 2, 1, 2, 3, 1, 0, 3, 0, 1, 1, 2, 3, 1, 1, 3, 2, 1, 2, 3, 0, 0))
prev_clust3c <- Rscc_clustering(c(2, 0, 3, 0, 2, 0, 2, 0, 2, 2, 1, 2, 3, 1, 0, 3, 0, 1, 1, 2, 3, 1, 1, 3, 2, 1, 2, 3, 0, 0))
prev_clust3d <- Rscc_clustering(c(2, 0, 3, 0, 2, 0, 2, 0, 2, 2, 1, 2, 3, 1, 0, 3, 0, 1, 1, 2, 3, 1, 1, 3, 2, 1, 2, 3, 0, 0))

test_that("greedy clustering functions cluster correctly without deep copy.", {
  expect_identical(top_down_greedy_clustering_internal(distance_object = test_distances2,
                                                       size_constraint = 2L,
                                                       batch_assign = TRUE,
                                                       existing_clustering = prev_clust3a,
                                                       deep_copy = FALSE),
                   replica_top_down_greedy_clustering(distance_object = test_distances2,
                                                      size_constraint = 2L,
                                                      batch_assign = TRUE,
                                                      existing_clustering = prev_clust3ref))
  expect_identical(top_down_greedy_clustering_internal(distance_object = test_distances2,
                                                       size_constraint = 3L,
                                                       batch_assign = TRUE,
                                                       existing_clustering = prev_clust3b,
                                                       deep_copy = FALSE),
                   replica_top_down_greedy_clustering(distance_object = test_distances2,
                                                      size_constraint = 3L,
                                                      batch_assign = TRUE,
                                                      existing_clustering = prev_clust3ref))
  expect_identical(top_down_greedy_clustering_internal(distance_object = test_distances2,
                                                       size_constraint = 2L,
                                                       batch_assign = FALSE,
                                                       existing_clustering = prev_clust3c,
                                                       deep_copy = FALSE),
                   replica_top_down_greedy_clustering(distance_object = test_distances2,
                                                      size_constraint = 2L,
                                                      batch_assign = FALSE,
                                                      existing_clustering = prev_clust3ref))
  expect_identical(top_down_greedy_clustering_internal(distance_object = test_distances2,
                                                       size_constraint = 3L,
                                                       batch_assign = FALSE,
                                                       existing_clustering = prev_clust3d,
                                                       deep_copy = FALSE),
                   replica_top_down_greedy_clustering(distance_object = test_distances2,
                                                      size_constraint = 3L,
                                                      batch_assign = FALSE,
                                                      existing_clustering = prev_clust3ref))
})

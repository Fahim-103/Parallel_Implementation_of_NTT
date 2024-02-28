#ifndef __TWIDDLE_H__
#define __TWIDDLE_H__


int32 chess_storage(%chess_alignof(int32)) stage0_tw[768] = { 1, 2352, 1854, 10302, 8685, 2802, 3400, 8950, 11632, 3150, 10822, 2825, 8340, 2436, 2798, 6281, 1534, 7291, 5277, 11903, 1514, 9407, 5064, 2487, 12149, 2523, 10798, 7822, 711, 968, 3271, 478, 5957, 1404, 8756, 10037, 12144, 3052, 1528, 5468, 6442, 11536, 10849, 4884, 9242, 10232, 3802, 8201, 7311, 3161, 12116, 10930, 11061, 11948, 9042, 6814, 1672, 64, 3060, 8055, 8011, 2835, 7282, 8687, 7506, 7108, 4976, 4424, 8754, 5333, 8436, 7026, 8736, 12153, 11931, 5925, 12163, 10873, 12176, 4582, 11700, 3329, 1715, 2888, 9048, 8637, 507, 431, 6014, 289, 3833, 7379, 3340, 3009, 10993, 11769, 5860, 6751, 964, 6152, 5351, 1616, 3531, 9837, 8726, 922, 5680, 1217, 11336, 7431, 2754, 1105, 5981, 8696, 4096, 11505, 11671, 8855, 9394, 11355, 2963, 1113, 219, 11239, 489, 7251, 9509, 11477, 7260, 6099, 3585, 1666, 10530, 4225, 7688, 5057, 10601, 11460, 4143, 11448, 497, 1489, 12052, 7870, 3006, 3937, 6207, 11821, 5274, 4847, 8241, 3079, 3587, 6370, 1949, 251, 480, 10661, 5112, 4782, 2829, 5459, 9852, 7139, 4154, 453, 8602, 4210, 9275, 1825, 3539, 4075, 11269, 9604, 1426, 11344, 1669, 5297, 9787, 1727, 6534, 6718, 9371, 6415, 9477, 9947, 9377, 8238, 8312, 10314, 42, 472, 4134, 2569, 8389, 7083, 7621, 7230, 9273, 9410, 12120, 8049, 6188, 4000, 6915, 5733, 2983, 11286, 432, 8366, 2143, 1846, 3775, 6142, 6409, 7654, 11112, 9010, 5284, 3789, 2203, 7787, 4414, 9812, 11371, 3728, 6199, 5294, 2731, 8454, 206, 5241, 965, 8504, 7205, 11918, 12216, 350, 12126, 9872, 5023, 4367, 9869, 10256, 11094, 3541, 8779, 2688, 5630, 6507, 4659, 8469, 10908, 8473, 8027, 3600, 1, 1854, 8685, 3400, 11632, 10822, 8340, 2798, 1534, 5277, 1514, 5064, 12149, 10798, 711, 3271, 5957, 8756, 12144, 1528, 6442, 10849, 9242, 3802, 7311, 12116, 11061, 9042, 1672, 3060, 8011, 7282, 7506, 4976, 8754, 8436, 8736, 11931, 12163, 12176, 11700, 1715, 9048, 507, 6014, 3833, 3340, 10993, 5860, 964, 5351, 3531, 8726, 5680, 11336, 2754, 5981, 4096, 11671, 9394, 2963, 219, 489, 9509, 7260, 3585, 10530, 7688, 10601, 4143, 497, 12052, 3006, 6207, 5274, 8241, 3587, 1949, 480, 5112, 2829, 9852, 4154, 8602, 9275, 3539, 11269, 1426, 1669, 9787, 6534, 9371, 9477, 9377, 8312, 42, 4134, 8389, 7621, 9273, 12120, 6188, 6915, 2983, 432, 2143, 3775, 6409, 11112, 5284, 2203, 4414, 11371, 6199, 2731, 206, 965, 7205, 12216, 12126, 5023, 9869, 11094, 8779, 5630, 4659, 10908, 8027, 79, 11287, 10220, 10531, 9542, 6997, 7543, 12129, 10585, 11346, 9005, 6808, 1229, 5101, 7013, 340, 3621, 3540, 834, 10111, 5069, 9130, 5067, 5422, 12275, 10911, 1300, 1556, 9198, 8249, 6130, 9984, 3102, 12145, 3382, 2838, 1960, 8585, 2335, 3362, 2625, 306, 2030, 3186, 8124, 7871, 5791, 8217, 8247, 2422, 4903, 8591, 1170, 6316, 10736, 8653, 5517, 4070, 334, 4786, 586, 5012, 1764, 1582, 8246, 568, 8507, 5191, 1827, 7783, 2396, 5855, 3983, 11082, 11109, 12011, 726, 6503, 1053, 10600, 2289, 4101, 8652, 3663, 7674, 9223, 5443, 2053, 8961, 11255, 48, 2969, 11343, 3443, 5331, 3318, 7072, 11414, 12187, 7516, 11227, 9581, 5569, 2166, 9550, 9540, 3289, 2462, 5329, 11899, 1991, 4614, 1212, 10450, 6836, 3985, 2501, 3901, 6522, 11701, 3569, 5444, 3907, 5357, 2366, 11680, 1502, 7394, 1, 10302, 3400, 3150, 8340, 6281, 5277, 9407, 12149, 7822, 3271, 1404, 12144, 5468, 10849, 10232, 7311, 10930, 9042, 64, 8011, 8687, 4976, 5333, 8736, 5925, 12176, 3329, 9048, 431, 3833, 3009, 5860, 6152, 3531, 922, 11336, 1105, 4096, 8855, 2963, 11239, 9509, 6099, 10530, 5057, 4143, 1489, 3006, 11821, 8241, 6370, 480, 4782, 9852, 453, 9275, 4075, 1426, 5297, 6534, 6415, 9377, 10314, 4134, 7083, 9273, 8049, 6915, 11286, 2143, 6142, 11112, 3789, 4414, 3728, 2731, 5241, 7205, 350, 5023, 10256, 8779, 6507, 10908, 3600, 11287, 156, 9542, 1973, 12129, 10695, 9005, 12138, 5101, 2738, 3621, 6427, 10111, 1958, 5067, 8851, 10911, 9928, 9198, 9606, 9984, 8527, 3382, 2049, 8585, 11026, 2625, 6950, 3186, 10542, 5791, 8076, 2422, 4774, 1170, 10120, 8653, 11089, 334, 12237, 5012, 7535, 8246, 8724, 5191, 8243, 2396, 7280, 11082, 1954, 726, 7540, 10600, 1146, 8652, 787, 9223, 9087, 8961, 1254, 2969, 11606, 5331, 421, 11414, 5876, 11227, 8775, 2166, 9597, 3289, 2505, 11899, 723, 1212, 400, 3985, 8210, 6522, 5681, 5444, 9381, 2366, 5445, 7394, 5766, 8595, 3445, 12047, 1583, 563, 11907, 9405, 3834, 1022, 9260, 9302, 11871, 7203, 4324, 10512, 3956, 4388, 6234, 354, 9364, 11567, 9090, 3000, 11454, 130, 12048, 11885, 3963, 2768, 5456, 10115, 6299, 6378, 9162, 7404, 10474, 5728, 10367, 9424, 2948, 4177, 7665, 8005, 8320, 9154, 11011, 7852, 5106, 5092, 8332, 9888, 2655, 8785, 6874, 6730, 10211, 12171, 975, 4337, 9259, 11289, 8471, 4053, 8273, 4231, 10968, 7270, 6374, 4821, 6093, 10163, 9235, 9821, 605, 2187, 4737, 955, 7210, 2704, 9734, 1428, 1323, 1045, 426	} ;



int32 chess_storage(%chess_alignof(int32)) stage1_tw[192]  = {  1, 8685, 11632, 8340, 1534, 1514, 12149, 711, 5957, 12144, 6442, 9242, 7311, 11061, 1672, 8011, 7506, 8754, 8736, 12163, 11700, 9048, 6014, 3340, 5860, 5351, 8726, 11336, 5981, 11671, 2963, 489, 7260, 10530, 10601, 497, 3006, 5274, 3587, 480, 2829, 4154, 9275, 11269, 1669, 6534, 9477, 8312, 4134, 7621, 12120, 6915, 432, 3775, 11112, 2203, 11371, 2731, 965, 12216, 5023, 11094, 5630, 10908, 1, 11632, 1534, 12149, 5957, 6442, 7311, 1672, 7506, 8736, 11700, 6014, 5860, 8726, 5981, 2963, 7260, 10601, 3006, 3587, 2829, 9275, 1669, 9477, 4134, 12120, 432, 11112, 11371, 965, 5023, 5630, 79, 9542, 10585, 1229, 3621, 5069, 12275, 9198, 3102, 1960, 2625, 8124, 8247, 1170, 5517, 586, 8246, 1827, 3983, 726, 2289, 7674, 8961, 11343, 7072, 11227, 9550, 5329, 1212, 2501, 3569, 2366, 1, 8340, 12149, 12144, 7311, 8011, 8736, 9048, 5860, 11336, 2963, 10530, 3006, 480, 9275, 6534, 4134, 6915, 11112, 2731, 5023, 10908, 9542, 9005, 3621, 5067, 9198, 3382, 2625, 5791, 1170, 334, 8246, 2396, 726, 8652, 8961, 5331, 11227, 3289, 1212, 6522, 2366, 8595, 563, 1022, 7203, 4388, 11567, 130, 2768, 6378, 5728, 4177, 9154, 5092, 8785, 12171, 11289, 4231, 4821, 9821, 955, 1428};

int32 chess_storage(%chess_alignof(int32)) stage2_tw[48]  = { 1, 1534, 5957, 7311, 7506, 11700, 5860, 5981, 7260, 3006, 2829, 1669, 4134, 432, 11371, 5023, 1, 5957, 7506, 5860, 7260, 2829, 4134, 11371, 79, 3621, 3102, 8247, 8246, 2289, 7072, 1212, 1, 7311, 5860, 3006, 4134, 5023, 3621, 2625, 8246, 8961, 1212, 563, 11567, 5728, 8785, 4821};

#endif
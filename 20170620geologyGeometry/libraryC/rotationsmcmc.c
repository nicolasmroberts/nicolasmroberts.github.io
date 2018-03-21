


/*
Copyright 2016 Joshua R. Davis

Licensed under the Apache License, Version 2.0 (the "License"); you may not 
use this file except in compliance with the License. You may obtain a copy of 
the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software 
distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
License for the specific language governing permissions and limitations under 
the License.
*/



/// COMPUTING THE POSTERIOR ///

// Lookup table for the function wtndJ below.
// Derived from numerical integration in Mathematica? I don't remember.
double wtndJData[301] = {
2.44948974278319, 2.4494897427831823, 2.449489742784173,
2.4494897427846993, 2.449489742784287, 2.4494897427834514, 
2.4494897427831623, 2.4494897427825006, 2.449489742781914, 
2.449489742780855, 2.449489742778968, 2.449489742775652, 
2.4494897427698765, 2.4494897427599396, 2.449489742743034, 
2.4494897427146083, 2.4494897428689617, 2.4494897425893503, 
2.4494897424623963, 2.4494897422576303, 2.449489741930722, 
2.449489741413934, 2.44948974060486, 2.449489739350127, 
2.4494897374223172, 2.449489734487214, 2.449489730058264, 
2.4494897234334485, 2.4494897136088, 2.4494896991609276, 
2.449489678088889, 2.4494896476032255, 2.449489603846946, 
2.4494895415296556, 2.449489453451811, 2.4494893298912523, 
2.4494891578185567, 2.449488919901518, 2.4494885932520574, 
2.449488148442007, 2.44948754465978, 2.449486733132628, 
2.4494856484071246, 2.449484207724826, 2.4494823061977264, 
2.449479811740776, 2.449476559063939, 2.4494723425998624, 
2.449466908236953, 2.449459943723486, 2.4494510676061334, 
2.4494398165668754, 2.449425631025655, 2.4494078388831166, 
2.4493856372882967, 2.449358072330842, 2.4493240165760533, 
2.4492821443841755, 2.449230904982693, 2.449168493291793, 
2.4490928185385115, 2.4490014707337107, 2.4488916851277054, 
2.4487602940784576, 2.44860371908509, 2.4484179241443305, 
2.4481983139798658, 2.4479397476402633, 2.447636509554559, 
2.447282231475882, 2.44686986180605, 2.4463916242962958, 
2.4458389781855363, 2.4452025802943567, 2.4444722495730278, 
2.4436369346066256, 2.4426846845752572, 2.4416026241230906, 
2.440376932579068, 2.438992827909754, 2.437434555752843, 
2.435685383811369, 2.433727601831453, 2.4315425273163846, 
2.42911051705808, 2.4264109844922666, 2.423422422807823, 
2.4201224336648854, 2.416487761301982, 2.412494331740721, 
2.408117296728728, 2.4033310819985174, 2.398109439362711, 
2.3924255021152376, 2.386251843164524, 2.3795605352884905, 
2.3723232128046026, 2.3645111343288034, 2.3560952452410744, 
2.347046240022701, 2.3373346232463605, 2.3269307687919216, 
2.315804970950179, 2.3039275157240504, 2.291268708827356, 
2.2777989350797347, 2.263488698267841, 2.248308660355831, 
2.2322296779424082, 2.2152228373346294, 2.1972594903220877, 
2.1783112933957147, 2.1583526993775606, 2.1373550636179663, 
2.115287258033935, 2.0921250173874792, 2.0678419949512286, 
2.0424122299414678, 2.0158109734156517, 1.9880139122309621, 
1.9589971552333167, 1.9287389060569982, 1.8972186040825614, 
1.8644158312555041, 1.830319825635167, 1.794912112092698, 
1.7581852930892967, 1.720135004839923, 1.6807628635344731, 
1.6400779226659448, 1.598098270030534, 1.5548529517645397, 
1.5103843903617122, 1.4647509793704603, 1.418030208635617, 
1.3703222415447018, 1.3217540321823484, 1.27248401999549, 
1.222707473273619, 1.1726625115929947, 1.1226368055037286, 
1.0729748590083688, 1.0240856070942153, 0.9764497509069772, 
0.9306257432572937, 0.8872525675649494, 0.8470464305582064, 
0.8107873785143679, 0.7792911950334243, 0.7533627101941639, 
0.733730089439142, 0.7209663764263263, 0.7154129006649309, 
0.7171246308831, 0.7258549331913795, 0.7410854293631779, 
0.7620914834293802, 0.7880238992694862, 0.8179871178654496, 
0.8511010447608743, 0.8865423422113804, 0.9235670146630912, 
0.9615190072215152, 0.9998297156648092, 1.038012374329191, 
1.0756540618877986, 1.1124070147290113, 1.1479801835206789, 
1.1821314857368188, 1.2146609228894942, 1.2454045779648597, 
1.2742294318024832, 1.301028930770894, 1.3257191516594087, 
1.3482355765809972, 1.368530298900488, 1.3865696535785663, 
1.402332194952716, 1.4158069770714774, 1.4269920964141871, 
1.4358934630736078, 1.4425237715491974, 1.4469016463902094, 
1.4490509412217318, 1.4490001723385857, 1.4467820701885377, 
1.4424332337971377, 1.435993874516497, 1.4275076367864827, 
1.4170214841736877, 1.4045856400058263, 1.3902535643191796, 
1.3740820134516705, 1.3561310045491937, 1.3364639566228322, 
1.31514771875909, 1.2922526454796663, 1.2678526550200737, 
1.2420252708106025, 1.2148516388394974, 1.1864165140696357, 
1.1568082096627306, 1.126118503450709, 1.0944424968909126, 
1.061878422653314, 1.0285273980125247, 0.9944931223475323, 
0.9598815182742396, 0.9248003172309595, 0.8893585916796297, 
0.853666237445236, 0.8178334110585664, 0.7819699282550034, 
0.7461846309755792, 0.7105847312849384, 0.6752751414565501, 
0.6403578003424288, 0.6059310062742533, 0.5720887673462763, 
0.5389201795533427, 0.5065088431430407, 0.474932326877238, 
0.44426168911634056, 0.4145610636068023, 0.38588731663890635, 
0.35828978088002306, 0.33181006971856514, 0.30648197441895425, 
0.28233144483085837, 0.25937665285802364, 0.23762813641236552, 
0.21708902019169043, 0.19775530835469898, 0.17961624304878418, 
0.16265472179354074, 0.14684776594893292, 0.13216703190895174, 
0.11857935626333374, 0.10604732595800898, 0.09452986445583975, 
0.08398282504071038, 0.07435958271216718, 0.06561161656520534, 
0.05768907512875261, 0.050541317824345684, 0.04411742648748594, 
0.038366681746076915, 0.03323899995291447, 0.028685327301998357,
0.024657988700572026, 0.021110989899354976, 0.018000272287447036,
0.015283920529780786, 0.012922324602027852, 0.010878296698172276, 
0.009117147192826492, 0.007606721201618443, 0.00631739968264778,
0.0052220686621212925, 0.004296060502396209, 0.003517071240393667,
0.0028650580558537927, 0.002322120859036193, 0.0018723718357744272, 
0.0015017965992959843, 0.001198110280431527, 0.000950611617244988,
0.0007500377253620529, 0.0005884218758423756, 0.0004589562313569278, 
0.0003558611111911294, 0.0002742619632404579, 0.00021007498686541502,
0.00015990194639747544, 0.00012093430792278923, 0.00009086703646081148, 
0.0000678215511915828, 0.00005027781163383719, 0.00003701456832687232,
0.000027058139410080313, 0.000019637734826775094, 0.000014148007432188404, 
0.000010116696230284038, 0.000007179276116811184, 0.0000050552290196200245,
0.0000035316648646876536, 0.0000024484813614844838, 0.0000016843116029229026, 
0.0000011492955894647053, 0.0000007800275572371351, 0.0000005277188862588222,
0.0000003565641869772526, 0.00000025685709711949437, 0.00000018311788389049094, 
0.00000014536004644075125, 0.00000012230311877112597, 0.00000011633584104021684,
0.00000009404267188322051, 0.00000009310693166331005, 0.00000007287508571373123, 0.0};

// Linearly interpolates the array wtndJData, which contains values for
// the Jeffreys prior J for the wrapped trivariate normal distribution 
// concentration parameter eta = -log kappa, from -1.0 to 2.0, sampled 
// at 0.01 intervals. See Qiu et al. (2014).
double wtndJ(double eta) {
  double i, t, lo, hi;
  int floori;
  i = (eta - -1.0) / 0.01;
  if (i <= 0.0)
    return wtndJData[0];
  else if (i >= 300.0) 
    return wtndJData[300];
  else {
    floori = (int)floor(i);
    t = i - floori;
    lo = wtndJData[floori];
    hi = wtndJData[floori + 1];
    return lo + (hi - lo) * t;
  }
}

// This is the (Lebesgue) angular density for the wrapped trivariate normal
// distribution, denoted g_wTND in Qiu et al. (2014).
double wtndAngularDensity(double r, double kappa, unsigned int terms) {
  int m;
  double sq, g = 0.0;
  for (m = -(int)terms; m <= (int)terms; m += 1) {
    sq = 2.0 * m * M_PI - r;
    sq = sq * sq;
    g += sq * exp(-kappa * kappa * sq * 0.5);
  }
  return g * kappa * kappa * kappa / sqrt(2.0 * M_PI);
}

double wtndLikelihood(unsigned int n, double *data, double *s, double kappa, unsigned int terms) {
  double sT[9], sTR[9], alpha, tr, like = 1.0;
  int i;
  matrixTranspose(s, sT);
  for (i = 0; i < n; i += 1) {
    // alpha = rotationDistance(s, datum).
    matrixMatrixMultiply(sT, &(data[9 * i]), sTR);
    tr = trace(sTR);
    alpha = arcCos((tr - 1.0) * 0.5);
    like *= wtndAngularDensity(alpha, kappa, terms) / (3.0 - tr);
  }
  return like;
}

/*
double wtndPosterior(unsigned int n, double *data, double *s, double eta, unsigned int terms) {
  return wtndLikelihood(n, data, s, exp(-eta), terms) * wtndJ(eta);
}
*/

// Returns 1 if challenger won, and 0 if incumbent won.
// Input s can safely alias output out.
int wtndNewS(unsigned int n, double *data, double *s, double kappa, unsigned int terms, double nu, double *out) {
  double numer, denom, r, sNew[9];
  rotationWrappedTrivariateNormal(s, exp(-nu), sNew);
  // The J factor cancels in the ratio, so just use the likelihood part.
  numer = wtndLikelihood(n, data, sNew, kappa, terms);
  denom = wtndLikelihood(n, data, s, kappa, terms);
  r = numer / denom;
  if (r >= 1.0 || doubleUniform(0.0, 1.0) <= r) {
    matrixCopy(sNew, out);
    return 1;
  } else {
    matrixCopy(s, out);
    return 0;
  }
}

// Returns 1 if challenger won, and 0 if incumbent won.
int wtndNewEta(unsigned int n, double *data, double *s, double eta, unsigned int terms,
    double logGamma, double *out) {
  double numer, denom, r, etaNew;
  etaNew = doubleNormal(eta, exp(logGamma));
  numer = wtndLikelihood(n, data, s, exp(-etaNew), terms) * wtndJ(etaNew);
  denom = wtndLikelihood(n, data, s, exp(-eta), terms) * wtndJ(eta);
  r = numer / denom;
  if (r >= 1.0 || doubleUniform(0.0, 1.0) <= r) {
    *out = etaNew;
    return 1;
  } else {
    *out = eta;
    return 0;
  }
}



/// TUNING ///

// Here are the ns and etas that are used as benchmarks.
double wtndNs[5] = {10.0, 30.0, 100.0, 300.0, 1000.0};
double wtndEtas[6] = {-4.0, -3.0, -2.0, -1.0, 0.0, 1.0};

// Each combination of n and eta produces a variance. Except that some are missing, coded 100.0.
double wtndVars10[6] = {0.00041, 0.00342, 0.02454, 0.18789, 100.0, 100.0};
double wtndVars30[6] = {0.00049, 0.00361, 0.02683, 0.19832, 1.36068, 100.0};
double wtndVars100[6] = {100.0, 0.00366, 0.0274, 0.20171, 1.45318, 1.90361};
double wtndVars300[6] = {100.0, 100.0, 100.0, 0.2017, 1.44678, 1.9271};
double wtndVars1000[6] = {100.0, 100.0, 100.0, 0.20264, 100.0, 100.0};

// Each combination of n and eta produces a tuned nu, except for the ones where variance == 100.0.
double wtndNus10[6] = {-4.95, -4.1, -2.75, -1.95, 100.0, 100.0};
double wtndNus30[6] = {-5.45, -4.45, -3.4, -2.45, -1.3, 100.0};
double wtndNus100[6] = {100.0, -5.05, -4.1, -3.05, -1.75, -3.15};
double wtndNus300[6] = {100.0, 100.0, 100.0, -3.6, -2.35, -3.5};
double wtndNus1000[6] = {100.0, 100.0, 100.0, -4.2, 100.0, 100.0};

// Each combination of n and eta produces a tuned gamma, roughly based on Qiu et al. (2014, Table 2).
double wtndGammas10[6] = {0.4, 0.4, 0.4, 0.4, 0.7, 0.6};
double wtndGammas30[6] = {0.23, 0.23, 0.23, 0.23, 0.5, 0.5};
double wtndGammas100[6] = {0.13, 0.13, 0.13, 0.13, 0.13, 0.3};
double wtndGammas300[6] = {0.07, 0.07, 0.07, 0.07, 0.08, 0.2};
double wtndGammas1000[6] = {0.04, 0.04, 0.04, 0.04, 0.05, 0.1};

void wtndLookup(unsigned int n, double *xs, double missing, double x,
    unsigned int *a, unsigned int *b, double *t) {
  int first, last, i;
  // The valid xs are those between xs[first] and xs[last], inclusive.
  first = 0;
  while (xs[first] == missing)
    first += 1;
  last = n - 1;
  while (xs[last] == missing)
    last -= 1;
  // Clamp if necessary.
  if (x < xs[first]) {
    *a = first;
    *b = first;
    *t = 0.0;
  } else if (xs[last] < x) {
    *a = last;
    *b = last;
    *t = 0.0;
  } else
    // Typical case.
    for (i = first; i < last; i += 1)
      if (xs[i] <= x && x <= xs[i + 1]) {
        *a = i;
        *b = i + 1;
        *t = (x - xs[i]) / (xs[i + 1] - xs[i]);
      }
}
  
// Guess eta and nu based on my experiments and gamma based on Qiu et al. (2014, Table 2).
void wtndSeed(unsigned int n, double variance, double *eta, double *nu, double *gamma) {
  // Quadrangulate based on n and variance.
  unsigned int a, b, aa, ab, ba, bb;
  double t, at, bt, aeta, beta, anu, bnu, agamma, bgamma;
  double *varss[5] = {wtndVars10, wtndVars30, wtndVars100, wtndVars300, wtndVars1000};
  wtndLookup(5, wtndNs, 100.0, (double)n, &a, &b, &t);
  wtndLookup(6, varss[a], 100.0, variance, &aa, &ab, &at);
  wtndLookup(6, varss[b], 100.0, variance, &ba, &bb, &bt);
  // We have four probably-redundant guesses for eta. Linearly interpolate them.
  aeta = wtndEtas[aa] + at * (wtndEtas[ab] - wtndEtas[aa]);
  beta = wtndEtas[ba] + bt * (wtndEtas[bb] - wtndEtas[ba]);
  *eta = aeta + t * (beta - aeta);
  // Linearly interpolate nu.
  double *nuss[5] = {wtndNus10, wtndNus30, wtndNus100, wtndNus300, wtndNus1000};
  anu = nuss[a][aa] + at * (nuss[a][ab] - nuss[a][aa]);
  bnu = nuss[b][ba] + bt * (nuss[b][bb] - nuss[b][ba]);
  *nu = anu + t * (bnu - anu);
  // Linearly interpolate gamma.
  double *gammass[5] = {wtndGammas10, wtndGammas30, wtndGammas100, wtndGammas300, wtndGammas1000};
  agamma = gammass[a][aa] + at * (gammass[a][ab] - gammass[a][aa]);
  bgamma = gammass[b][ba] + bt * (gammass[b][bb] - gammass[b][ba]);
  *gamma = agamma + t * (bgamma - agamma);
  fprintf(stderr, "wtndSeed: n=%d, var=%f, a=%d, b=%d, aa=%d, ab=%d, ba=%d, bb=%d, eta=%f, nu=%f, gamma=%f\n",
    n, variance, a, b, aa, ab, ba, bb, *eta, *nu, *gamma);
}

// Tunes the parameter value to keep the acceptance rate between 0.30 and 0.40.
// The old value produced the old rate. The new value produced the rate countNew / tuning.
// General principle: If the acceptance rate is too high, then we are "playing
// it too safe", and we should go for larger steps. If the acceptance rate is 
// too low, then we should go for smaller steps. But this doesn't always pan out.
void wtndTuning(double *rateOld, double *paramOld, int tuning, int *countNew, double *paramNew) {
  double param, rateNew;
  rateNew = *countNew / (double)tuning;
  if (0.3 <= rateNew && rateNew <= 0.4)
    param = *paramNew;
  else if (rateNew < 0.3)
    param = *paramNew - doubleUniform(0.05, 0.15);
  else
    param = *paramNew + doubleUniform(0.05, 0.15);
  //fprintf(stderr, "rateOld %f, paramOld %f, rateNew %f, paramNew %f, param %f\n",
  //  *rateOld, *paramOld, rateNew, *paramNew, param);
  *rateOld = rateNew;
  *paramOld = *paramNew;
  *countNew = 0;
  *paramNew = param;
}



/// MARKOV CHAIN MONTE CARLO ///

// This function is used in experiments about acceptance rates.
// It doesn't alter any parameters except nuRate and gammaRate.
void mcmcRotationRate(unsigned int n, double *data, unsigned int tuning,
    double *s, double eta, double nu, double logGamma, double *nuRate, double *gammaRate) {
  int i, terms = 10, nuCount = 0, gammaCount = 0;
  double ss[9];
  matrixCopy(s, ss);
  for (i = 0; i < tuning; i += 1) {
    nuCount += wtndNewS(n, data, ss, exp(-eta), terms, nu, ss);
    gammaCount += wtndNewEta(n, data, ss, eta, terms, logGamma, &eta);
  }
  *nuRate = nuCount / (double)tuning;
  *gammaRate = gammaCount / (double)tuning;
}

void mcmcRotationSeeding(unsigned int n, double *data, unsigned int tuning,
    double *sSeed, double *etaSeed, double *nuSeed, double *gammaSeed) {
  // Guess at reasonable seed values for (S, eta) and (nu, gamma).
  double nuGuess, gammaGuess, nu, gamma, logGamma, error, variance;
  unsigned int used;
  rotationMeanMany(n, data, 10, 0.0000001, 1000, sSeed, &error, &used, &variance);
  wtndSeed(n, variance, etaSeed, &nuGuess, &gammaGuess);
  // This starting value for best error-squared is impossibly bad, so it will get overwritten.
  double nuRates[5], gammaRates[5];
  double nuRate, gammaRate, errSq, bestErrSq = 4.0;
  int i, j, k, progress = 0;
  // Try nus between nuSeed - 1.5 and nuSeed + 1.5.
  nu = nuGuess - 2.0;
  for (i = -3; i <= 3; i += 1) {
    nu += 0.5;
    // Try gammas between gammaSeed / 4 and gammaSeed * 4.
    gamma = gammaGuess * 0.125;
    for (j = -2; j <= 2; j += 1) {
      fprintf(stderr, "%f ", progress / 35.0);
      progress += 1;
      gamma *= 2.0;
      logGamma = log(gamma);
      for (k = 0; k < 5; k += 1)
        mcmcRotationRate(n, data, tuning, sSeed, *etaSeed, nu, logGamma, &(nuRates[k]), &(gammaRates[k]));
      nuRate = doubleMedian(5, nuRates);
      gammaRate = doubleMedian(5, gammaRates);
      // If these rates are closest to 0.35 of any seen yet, then remember them.
      errSq = (nuRate - 0.35) * (nuRate - 0.35) + (gammaRate - 0.35) * (gammaRate - 0.35);
      if (errSq < bestErrSq) {
        bestErrSq = errSq;
        *nuSeed = nu;
        *gammaSeed = gamma;
      }
    }
  }
  fprintf(stderr, "\n");
}

// data is an array of 9 * n doubles, to store the data rotations in column-major order.
// terms controls certain asymptotic expansions; 10 seems like a good value.
// tuning is the number of MCMC iterations between tunings; 1,000 or 10,000 seems good.
// burnin bounds the number of tunings undertaken in the burn-in phase;
// 100 or 1,000 seems good, for a total of about 10^6 MCMC iterations.
// collection is the number of tunings in the collection phase;
// 10^7 / tuning is common, for a total of 10^7 MCMC iterations.
// meta is an array of 8 doubles reporting nu, nuRate, logGamma, and gammaRate
// at the end of burn-in and the end of collection.
// sBuffer is an array of 9 * tuning * collection doubles, to store the Ss in column-major order.
// etaBuffer is an array of tuning * collection doubles, to store the etas.
// Before calling this function, ensure that initializeRandom has been called.
void mcmcRotationWrappedTrivariateNormal(unsigned int n, double *data,
    unsigned int terms, unsigned int tuning, unsigned int burnin, unsigned int collection,
    double *meta, double *sBuffer, double *etaBuffer) {
  double s[9], eta, nu, gamma, logGamma;
  fprintf(stderr, "mcmcRotationWrappedTrivariateNormal: beginning seeding at ");
  printTime(stderr);
  fprintf(stderr, "\n");
  mcmcRotationSeeding(n, data, tuning, s, &eta, &nu, &gamma);
  logGamma = log(gamma);
  fprintf(stderr, "mcmcRotationWrappedTrivariateNormal: beginning burn-in at ");
  printTime(stderr);
  fprintf(stderr, "\n");
  int i, j, stable = 0, nuCount = 0, gammaCount = 0;
  double nuRateOld = -1.0, nuOld, gammaRateOld = -1.0, logGammaOld;
  i = 0;
  while (i < burnin && stable < 3) {
    i += 1;
    for (j = 0; j < tuning; j += 1) {
      nuCount += wtndNewS(n, data, s, exp(-eta), terms, nu, s);
      gammaCount += wtndNewEta(n, data, s, eta, terms, logGamma, &eta);
    }
    wtndTuning(&nuRateOld, &nuOld, tuning, &nuCount, &nu);
    wtndTuning(&gammaRateOld, &logGammaOld, tuning, &gammaCount, &logGamma);
    if (nu == nuOld && logGamma == logGammaOld)
      stable += 1;
    else
      stable = 0;
    fprintf(stderr, "tuning %d, last nu-rate %f, last gamma-rate %f, nu %f, logGamma %f\n",
      i, nuRateOld, gammaRateOld, nu, logGamma);
  }
  // Store the first chunk of metadata.
  meta[0] = nu;
  meta[1] = nuRateOld;
  meta[2] = logGamma;
  meta[3] = gammaRateOld;
  if (0.3 <= nuRateOld && nuRateOld <= 0.4 && 0.3 <= gammaRateOld && gammaRateOld <= 0.4) {
    fprintf(stderr, "mcmcRotationWrappedTrivariateNormal: beginning collection at ");
    printTime(stderr);
    fprintf(stderr, "\n");
    for (i = 0; i < collection; i += 1) {
      for (j = 0; j < tuning; j += 1) {
        nuCount += wtndNewS(n, data, s, exp(-eta), terms, nu, s);
        gammaCount += wtndNewEta(n, data, s, eta, terms, logGamma, &eta);
        matrixCopy(s, &(sBuffer[9 * (i * tuning + j)]));
        etaBuffer[i * tuning + j] = eta;
      }
      wtndTuning(&nuRateOld, &nuOld, tuning, &nuCount, &nu);
      wtndTuning(&gammaRateOld, &logGammaOld, tuning, &gammaCount, &logGamma);
      fprintf(stderr, "tuning %d, last nu-rate %f, last gamma-rate %f, nu %f, logGamma %f\n",
        i, nuRateOld, gammaRateOld, nu, logGamma);
    }
    // Store the second chunk of metadata.
    meta[4] = nu;
    meta[5] = nuRateOld;
    meta[6] = logGamma;
    meta[7] = gammaRateOld;
  } else
    fprintf(stderr, "mcmcRotationWrappedTrivariateNormal: aborting after failed burn-in\n");
}

// This function is used in experiments about acceptance rates.
double mcmcRotationSurgical(unsigned int n, double *data, unsigned int tuning, double nuTest) {
  // Guess at reasonable values for (S, eta) and (nu, gamma).
  double s[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double error, variance, eta, nu, gamma, logGamma;
  unsigned int used;
  rotationMean(n, data, s, 0.0000001, 1000, s, &error, &used);
  variance = rotationVariance(n, data, s);
  wtndSeed(n, variance, &eta, &nu, &gamma);
  logGamma = log(gamma);
  // Prepare for loop.
  int i, terms = 10;
  int nuCount = 0, gammaCount = 0;
  // Run the MCMC but don't count yet.
  for (i = 0; i < tuning; i += 1) {
    nuCount += wtndNewS(n, data, s, exp(-eta), terms, nuTest, s);
    gammaCount += wtndNewEta(n, data, s, eta, terms, logGamma, &eta);
  }
  // Now count.
  nuCount = 0;
  gammaCount = 0;
  for (i = 0; i < tuning; i += 1) {
    nuCount += wtndNewS(n, data, s, exp(-eta), terms, nuTest, s);
    gammaCount += wtndNewEta(n, data, s, eta, terms, logGamma, &eta);
  }
  return nuCount / (double)tuning;
}



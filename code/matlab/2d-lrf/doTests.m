rep = 10;

Theta_tqr_mi = zeros(rep, 3);
Sigma_tqr_mi = zeros(3, 3, rep);
batchSize = zeros(rep, 1);

Theta_ls = zeros(rep, 3);
Sigma_ls = zeros(3, 3, rep);

Theta_ekf = zeros(rep, 3);
Sigma_ekf = zeros(3, 3, rep);

params;

for i = 1:rep
  % setting up the scene
  clear u;
  clear l;
  u = genStraightPath(5000, 0.1);
  l = genLandmarks(17, [0 40], [0 40]);
  [x_true y_true th_true r b v om t] = simulate(u, x0, l, Q, R, Theta, T);
  x_odom = odomInt([x_true(1), y_true(1), th_true(1)], [v, om], t);
  l_hat = initLandmarks(x_odom, Theta_hat, r, b);

  % TQR-MI optimization
  [x_est l_est Theta_est Sigma batchIdx] = ls_slam_calib_it(x_odom, l_hat, ...
    Theta_hat, u, r, b, t, Q, R, 50, 1e-9, 100, 0.4, 0.02);
  Theta_tqr_mi(i, :) = Theta_est';
  Sigma_tqr_mi(:, :, i) = Sigma;
  Theta_est
  Sigma
  batchSize(i) = length(batchIdx);

  % Non-linear least squares without regularization
  [x_est l_est Theta_est Sigma] = ls_slam_calib(x_odom, l_hat, Theta_hat, u, ...
    r, b, t, Q, R, 50, 1e-9);
  Theta_ls(i, :) = Theta_est';
  Sigma_ls(:, :, i) = Sigma;
  Theta_est
  Sigma

  % EKF
  [x_est, P] = ekf_slam_calib([x_true(1); y_true(1); th_true(1)], ...
    diag([1e-6; 1e-6; 1e-6]), Theta_hat, diag([0.0025; 0.0025; 0.0025]), u, ...
    r, b, t, Q, R);
  Theta_ekf(i, :) = x_est(end, end - 2:end);
  Sigma_ekf(:, :, i) = P(end - 2:end, end - 2:end, end);
  x_est(end, end - 2:end)
  P(end - 2:end, end - 2:end, end)

end

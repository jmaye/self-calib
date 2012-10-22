for i = 1:1
  % setting up the scene
  clear all;
  params;
  u = genStraightPath(5000,0.1);
  l = genLandmarks(17,[0 40],[0 40]);
  [x_true y_true th_true r b v om t] = simulate(u, x0, l, Q, R, Theta, T);
  x_odom = odomInt([x_true(1), y_true(1), th_true(1)], [v, om], t);
  l_hat = initLandmarks(x_odom, Theta_hat, r, b);

  % TQR-MI optimization
  [x_est l_est Theta_est Sigma batchIdx] = ls_slam_calib_it(x_odom, l_hat, ...
    Theta_hat, u, r, b, t, Q, R, 200, 1e-9, 100, 0.4, 0.018);
end

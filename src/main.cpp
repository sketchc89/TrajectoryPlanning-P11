#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "spline.h"
#include "json.hpp"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
constexpr double deg2rad(double x) { return x * pi() / 180; }
constexpr double rad2deg(double x) { return x * 180 / pi(); }
constexpr double mph2Mps(double x) { return x * 0.44704; }
constexpr double mps2Mph(double x) { return x * 2.23694; }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for(int i = 0; i < maps_x.size(); i++)
  {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if(dist < closestLen)
    {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y-y),(map_x-x));

  double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
    if (closestWaypoint == maps_x.size())
    {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
  {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
  {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}

bool laneChangeSafe(int current_lane, int target_lane, const vector<vector<double>>& vehicles) {

  cout << "\nCHECKING: Attempting change from lane " << current_lane << " to " << target_lane << "\n";
  if(current_lane == target_lane) {
    cout << "SAFE: Target lane same as current. No lane change needed\n";
    return true;
  }
  if(abs(current_lane - target_lane) >=2) {
    cout << "UNSAFE: Target lane more than one lane away. Not safe\n";
    return false;
  }

  double v_host = sqrt(vehicles[0][3]*vehicles[0][3] + vehicles[0][4]*vehicles[0][4]);
  double s_host = vehicles[0][5];
  bool ahead_of_host, behind_host, side_of_host;

  for(auto i=1; i<vehicles.size(); ++i) {
    double vx_actor = vehicles[i][3];
    double vy_actor = vehicles[i][4];
    double s_actor = vehicles[i][5];
    double d_actor = vehicles[i][6];
    //check if actor near target lane
    if(s_actor - s_host > 3) { 
      ahead_of_host = true; 
      behind_host = false; 
      side_of_host = false;
    } else if(s_host - s_actor > 3) {
      ahead_of_host = false;
      behind_host = true;
      side_of_host = false;
    } else {
      ahead_of_host = false;
      behind_host = false;
      side_of_host = true;
    }

    double clearance = v_actor; // TODO ignore vehicles outside clearance
    if(target_lane == 1 && d_actor > 4) { continue; }
    else if(target_lane == 2 && (d_actor < 4 || d_actor < 8)) { continue; }
    else if(target_lane == 3 && d_actor < 8) { continue; }
    double v_actor = sqrt(vx_actor*vx_actor + vy_actor*vy_actor);
    if (ahead_of_host && v_host > v_actor && s_host - s_actor < clearance) { //TODO fix
      cout << "UNSAFE: Actor vehicle ahead in target lane going too slow\n";
      cout << "\tV_host " << v_host << "\tV_actor " << v_actor << "\tS_host " << s_host << "\tS_actor " << s_actor << "\n";
      return false; 
    }
    if (side_of_host) { 
      cout << "Actor vehicle to side in target lane\n";
      return false; 
    }
    if (behind_host && s_host - s_actor < clearance && v_actor > v_host) { 
      cout << "Actor vehicle behind in target lane too close\n";
      cout << "\tV_host " << v_host << "\tV_actor " << v_actor << "\tS_host " << s_host << "\tS_actor " << s_actor << "\n";
      return false; }
  }

  cout << "SAFE\n";
  return true;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  double v_set = 0.0;
  h.onMessage([&v_set, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
        uWS::OpCode opCode) {
      // "42" at the start of the message means there's a websocket message event.
      // The 4 signifies a websocket message
      // The 2 signifies a websocket event
      //auto sdata = string(data).substr(0, length);
      //cout << sdata << endl;
      if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
      auto j = json::parse(s);

      string event = j[0].get<string>();

      if (event == "telemetry") {
      double x_car = j[1]["x"];
      double y_car = j[1]["y"];
      double s_car = j[1]["s"];
      double d_car = j[1]["d"];
      double yaw_car = j[1]["yaw"];
      double v_car = j[1]["speed"];

      // Previous path data given to the Planner
      auto previous_path_x = j[1]["previous_path_x"];
      auto previous_path_y = j[1]["previous_path_y"];
      // Previous path's end s and d values 
      double end_path_s = j[1]["end_path_s"];
      double end_path_d = j[1]["end_path_d"];

      // Sensor Fusion Data, a list of all other cars on the same side of the road.
      vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];
      vector<double> host_vehicle {-1, x_car, y_car, v_car*cos(yaw_car), v_car*sin(yaw_car), s_car, d_car};
      sensor_fusion.insert(sensor_fusion.begin(), host_vehicle);
      laneChangeSafe(1,1, sensor_fusion);
      laneChangeSafe(1,2, sensor_fusion);
      laneChangeSafe(1,3, sensor_fusion);
      laneChangeSafe(2,1, sensor_fusion);
      laneChangeSafe(2,2, sensor_fusion);
      laneChangeSafe(2,3, sensor_fusion);
      laneChangeSafe(3,1, sensor_fusion);
      laneChangeSafe(3,2, sensor_fusion);
      laneChangeSafe(3,3, sensor_fusion);
      /*for (auto i=0; i<sensor_fusion.size(); ++i) {
        cout << "ID\tX\tY\tVx\tVy\tS\tD\n";
        cout << sensor_fusion[i][0] 
          << "\t" << sensor_fusion[i][1]
          << "\t" << sensor_fusion[i][2]
          << "\t" << sensor_fusion[i][3]
          << "\t" << sensor_fusion[i][4]
          << "\t" << sensor_fusion[i][5]
          << "\t" << sensor_fusion[i][6] << "\n";
      }*/


      json msgJson;

      vector<double> next_x_vals, next_y_vals, x_anchors, y_anchors;
      int lane = 2;
      int max_path_size = 50;
      double v_ref = mph2Mps(49.0);
      double x_ref = x_car;
      double y_ref = y_car;
      double x_prev_ref, y_prev_ref;
      double yaw_ref = deg2rad(yaw_car);
      int prev_size = previous_path_x.size();

      double v_set;
      double tolerance = 0.4;
      if(v_set < v_ref - tolerance) { // speed up
        v_set += 0.4;
      } else  if(v_set > v_ref + tolerance) {
        v_set -= 0.4;
      } // else stay the same

      cout << "V_ref: " << v_ref << "\tX_ref: " << x_ref <<  "\tY_ref: " << y_ref <<  "\tYaw_ref: " << yaw_ref << "\tPrev_size: " << prev_size << "\n";

      if(prev_size <2) {
        x_prev_ref = x_ref - cos(yaw_ref);
        y_prev_ref = y_ref - sin(yaw_ref);
        cout << "no prev path, adding points x,y " << x_prev_ref << ", " << y_prev_ref << "\t" << x_ref << ", " << y_ref << "\n";
        x_anchors.push_back(x_prev_ref);
        x_anchors.push_back(x_ref);
        y_anchors.push_back(y_prev_ref);
        y_anchors.push_back(y_ref);
      } else {
        x_ref =      previous_path_x[prev_size-1];
        x_prev_ref = previous_path_x[prev_size-2];
        y_ref =      previous_path_y[prev_size-1];
        y_prev_ref = previous_path_y[prev_size-2];
        yaw_ref = atan2(y_ref - y_prev_ref, x_ref - x_prev_ref);

        x_anchors.push_back(x_prev_ref);
        x_anchors.push_back(x_ref);
        y_anchors.push_back(y_prev_ref);
        y_anchors.push_back(y_ref);
      }
      vector<double> next_wp0, next_wp1, next_wp2;
      next_wp0 = getXY(s_car+30, 4*(lane-0.5), map_waypoints_s, map_waypoints_x, map_waypoints_y);
      next_wp1 = getXY(s_car+60, 4*(lane-0.5), map_waypoints_s, map_waypoints_x, map_waypoints_y);
      next_wp2 = getXY(s_car+90, 4*(lane-0.5), map_waypoints_s, map_waypoints_x, map_waypoints_y);

      x_anchors.push_back(next_wp0[0]);
      x_anchors.push_back(next_wp1[0]);
      x_anchors.push_back(next_wp2[0]);
      y_anchors.push_back(next_wp0[1]);
      y_anchors.push_back(next_wp1[1]);
      y_anchors.push_back(next_wp2[1]);

      for(auto i=0; i< x_anchors.size(); ++i) {
        double x_shift = x_anchors[i] - x_ref;
        double y_shift = y_anchors[i] - y_ref;
        x_anchors[i] = x_shift*cos(-yaw_ref) - y_shift*sin(-yaw_ref);
        y_anchors[i] = x_shift*sin(-yaw_ref) + y_shift*cos(-yaw_ref);
        cout << "I: " << i << "\tX: " << x_anchors[i] << "\tY: " << y_anchors[i] << "\n";
      }
      tk::spline s;
      s.set_points(x_anchors, y_anchors);

      for(auto i=0; i<prev_size; ++i) {
        next_x_vals.push_back(previous_path_x[i]);
        next_y_vals.push_back(previous_path_y[i]);
      }

      double x_target = 30.0;
      double y_target = s(x_target);
      double dist_target = sqrt( x_target*x_target + y_target*y_target);
      double x_addon = 0;
      for(auto i=0; i<max_path_size - prev_size; ++i) {
        double N = dist_target/(0.02*v_set);
        double x_pt = x_addon + x_target/N;
        double y_pt = s(x_pt);
        x_addon = x_pt;
        double x_temp = x_pt;
        double y_temp = y_pt;
        //cout << "X_temp: " << x_temp << "\tY_temp: " << y_temp << "\n";

        x_pt = x_ref + x_temp*cos(yaw_ref) - y_temp*sin(yaw_ref);
        y_pt = y_ref + x_temp*sin(yaw_ref) + y_temp*cos(yaw_ref);
        //cout << "X_pt: " << x_pt << "\tY_pt: " << y_pt << "\n";
        next_x_vals.push_back(x_pt);
        next_y_vals.push_back(y_pt);
      }

      msgJson["next_x"] = next_x_vals;
      msgJson["next_y"] = next_y_vals;

      auto msg = "42[\"control\","+ msgJson.dump()+"]";

      //this_thread::sleep_for(chrono::milliseconds(1000));
      ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

      }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
      }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
        size_t, size_t) {
      const std::string s = "<h1>Hello world!</h1>";
      if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
      } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
      }
      });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
      std::cout << "Connected!!!" << std::endl;
      });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
        char *message, size_t length) {
      ws.close();
      std::cout << "Disconnected" << std::endl;
      });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

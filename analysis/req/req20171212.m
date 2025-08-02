% req20171212 ----------------------------------------------------
%  Status: active
%  Type: ProcessData
%  Final_Forms: NA
%  Description: Go through Eduardo's rats which contain pressure 
%               sensors to collect statistics.
%  Bugs: NA

Trial = MTATrial.validate('Ed01-20140707.cof.all'); % GOOD
Trial = MTATrial.validate('Ed01-20140709.cof.all'); % GOOD
Trial = MTATrial.validate('Ed01-20140717.cof.all'); % GOOD

Trial = MTATrial.validate('Ed03-20140624.cof.all'); % WEAK
Trial = MTATrial.validate('Ed03-20140625.cof.all'); % WEAK

Trial = MTATrial.validate('Ed05-20140528.cof.all'); % GOOD
Trial = MTATrial.validate('Ed05-20140529,ont.all'); % GOOD

Trial = MTASession.validate('Ed10-20140812.cof.all');
Trial = MTASession.validate('Ed10-20140814.cof.all');
Trial = MTASession.validate('Ed10-20140812.cof.all');
Trial = MTATrial.validate('Ed10-20140812.cof.all');
Ed10-20140814
Ed10-20140815
Ed10-20140816
Ed10-20140817

Ed11-20150807
Ed11-20150811
Ed11-20150813
Ed11-20151006
Ed11-20151016
Ed11-20151020
Ed11-20151021
Ed11-20151104
Ed11-20151108

Ed12-20150807
Ed12-20150810
Ed12-20150811
Ed12-20150813
Ed12-20151014
Ed12-20151019
Ed12-20151020
Ed12-20151021
Ed12-20151104
Ed12-20151108



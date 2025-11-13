function out = model
%
% jqm.m
%
% Model exported on Nov 11 2025, 23:56 by COMSOL 6.3.0.290.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath(['C:\Users\' native2unicode(hex2dec({'51' 'b7'}), 'unicode')  native2unicode(hex2dec({'51' '74'}), 'unicode') '\Desktop']);

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').physics.create('spf', 'LaminarFlow', 'geom1');
model.component('comp1').physics.create('pf', 'PhaseField', 'geom1');

model.component('comp1').multiphysics.create('tpf1', 'TwoPhaseFlowPhaseField', 2);
model.component('comp1').multiphysics('tpf1').set('Fluid_physics', 'spf');
model.component('comp1').multiphysics('tpf1').set('Mathematics_physics', 'pf');
model.component('comp1').multiphysics('tpf1').selection.all;

model.study.create('std1');
model.study('std1').create('phasei', 'PhaseInitialization');
model.study('std1').feature('phasei').set('ftplistmethod', 'manual');
model.study('std1').feature('phasei').set('solnum', 'auto');
model.study('std1').feature('phasei').set('notsolnum', 'auto');
model.study('std1').feature('phasei').set('outputmap', {});
model.study('std1').feature('phasei').set('ngenAUX', '1');
model.study('std1').feature('phasei').set('goalngenAUX', '1');
model.study('std1').feature('phasei').set('ngenAUX', '1');
model.study('std1').feature('phasei').set('goalngenAUX', '1');
model.study('std1').feature('phasei').setSolveFor('/physics/spf', true);
model.study('std1').feature('phasei').setSolveFor('/physics/pf', true);
model.study('std1').feature('phasei').setSolveFor('/multiphysics/tpf1', true);
model.study('std1').feature('phasei').setSolveFor('/physics/spf', false);
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('initstudy', 'std1');
model.study('std1').feature('time').set('notstudy', 'std1');
model.study('std1').feature('time').set('ftplistmethod', 'manual');
model.study('std1').feature('time').set('initialtime', '0');
model.study('std1').feature('time').set('useinitsol', 'on');
model.study('std1').feature('time').set('notsolmethod', 'sol');
model.study('std1').feature('time').set('outputmap', {});
model.study('std1').feature('time').setSolveFor('/physics/spf', true);
model.study('std1').feature('time').setSolveFor('/physics/pf', true);
model.study('std1').feature('time').setSolveFor('/multiphysics/tpf1', true);

model.component('comp1').geom('geom1').lengthUnit('Mm');
model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').run('r2');
model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').feature('r1').set('size', [1 5.35]);
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').feature('r2').set('size', [5 10]);
model.component('comp1').geom('geom1').feature('r2').set('base', 'center');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'0' '010-5.35/2+0.35'});
model.component('comp1').geom('geom1').run('r2');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'0' '10-5.35/2+0.35'});
model.component('comp1').geom('geom1').run('r2');
model.component('comp1').geom('geom1').feature('r3').set('size', [4 2.4]);
model.component('comp1').geom('geom1').feature('r3').set('pos', {'0' '10-5.35/2+0.35'});
model.component('comp1').geom('geom1').run('r3');
model.component('comp1').geom('geom1').feature('r3').set('base', 'center');
model.component('comp1').geom('geom1').run('r3');
model.component('comp1').geom('geom1').feature('r3').set('pos', {'0' '1.675-0.2'});
model.component('comp1').geom('geom1').run('r3');
model.component('comp1').geom('geom1').create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol1').set('table', [-2 0.275; -0.5 -2.675; 0.5 -2.675; 2 0.275; -2 0.275]);
model.component('comp1').geom('geom1').run('pol1');
model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').set('size', [40 18]);
model.component('comp1').geom('geom1').feature('r4').set('pos', {'0-20' '0-10'});
model.component('comp1').geom('geom1').run('r4');
model.component('comp1').geom('geom1').create('uni1', 'Union');
model.component('comp1').geom('geom1').feature('uni1').selection('input').set({'pol1' 'r1' 'r3'});
model.component('comp1').geom('geom1').run('uni1');
model.component('comp1').geom('geom1').feature('uni1').set('intbnd', false);
model.component('comp1').geom('geom1').run('uni1');
model.component('comp1').geom('geom1').feature.duplicate('r5', 'r4');
model.component('comp1').geom('geom1').feature('r5').set('size', [40 8]);
model.component('comp1').geom('geom1').run('r5');
model.component('comp1').geom('geom1').create('uni2', 'Union');
model.component('comp1').geom('geom1').feature('uni2').selection('input').set({'r2' 'r4'});
model.component('comp1').geom('geom1').run('uni2');
model.component('comp1').geom('geom1').run('r2');
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').run('r2');
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').run('r2');
model.component('comp1').geom('geom1').run('r3');
model.component('comp1').geom('geom1').run('r5');
model.component('comp1').geom('geom1').run('uni2');
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').create('mce1', 'MeshControlEdges');
model.component('comp1').geom('geom1').feature('mce1').selection('input').set('fin', [4 15 21]);
model.component('comp1').geom('geom1').run('mce1');

model.component('comp1').common.create('free1', 'DeformingDomain');
model.component('comp1').common('free1').selection.all;
model.component('comp1').common('free1').selection.set([1 4]);
model.component('comp1').common.create('pnmd1', 'PrescribedNormalMeshDisplacement');
model.component('comp1').common('pnmd1').selection.set([3 20]);
model.component('comp1').common.create('pres1', 'PrescribedDeformation');
model.component('comp1').common('pres1').selection.set([2 3 4]);
model.component('comp1').common('pres1').set('prescribedDeformation', {'10[mm]*sin(0.5*2*pi*t[1/s])' '0' '0'});

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an1', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an2', 'Analytic');
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup.create('NonlinearModel', 'NonlinearModel', 'Nonlinear model');
model.component('comp1').material('mat1').propertyGroup.create('idealGas', 'idealGas', 'Ideal gas');
model.component('comp1').material('mat1').propertyGroup('idealGas').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat1').label('Air');
model.component('comp1').material('mat1').set('family', 'air');
model.component('comp1').material('mat1').propertyGroup('def').label('Basic');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').label('Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').label('Piecewise 2');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').label('Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/R_const[K*mol/J]/T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', {'Pa' 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotaxis', {'off' 'on'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotfixedvalue', {'101325' '273.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '293.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').label('Piecewise 3');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').label('Analytic 2');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*R_const[K*mol/J]/0.02897*T)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', 'm/s');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotfixedvalue', {'273.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').label('Analytic 1');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(pA,T)*d(rho(pA,T),T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'pA' 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', {'Pa' 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotaxis', {'off' 'on'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotfixedvalue', {'101325' '273.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').label('Analytic 2a');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'muB');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '0.6*eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotfixedvalue', {'200'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T' '200' '1600'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '');
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '0.02897[kg/mol]');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(pA,T)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('pressure');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').label('Refractive index');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('NonlinearModel').label('Nonlinear model');
model.component('comp1').material('mat1').propertyGroup('NonlinearModel').set('BA', 'def.gamma-1');
model.component('comp1').material('mat1').propertyGroup('idealGas').label('Ideal gas');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').label('Piecewise 2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('Rs', 'R_const/Mn');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('molarmass', '0.02897[kg/mol]');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('pressure');
model.component('comp1').material('mat1').materialType('nonSolid');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').label([native2unicode(hex2dec({'9a' 'd8'}), 'unicode')  native2unicode(hex2dec({'7c' '98'}), 'unicode')  native2unicode(hex2dec({'5e' 'a6'}), 'unicode')  native2unicode(hex2dec({'6d' '41'}), 'unicode')  native2unicode(hex2dec({'4f' '53'}), 'unicode') ]);
model.component('comp1').material('mat2').selection.set([2]);
model.component('comp1').material('mat2').propertyGroup('def').set('density', {'1750'});
model.component('comp1').material('mat2').propertyGroup('def').set('dynamicviscosity', {'2.35'});

model.component('comp1').physics('spf').prop('PhysicalModelProperty').set('IncludeGravity', true);
model.component('comp1').physics('spf').create('open1', 'OpenBoundary', 1);
model.component('comp1').physics('spf').feature('open1').selection.set([1 3 20 21]);
model.component('comp1').physics('spf').create('iwbc1', 'InteriorWallBC', 1);
model.component('comp1').physics('spf').feature('iwbc1').selection.set([4 5 9 10 14 16 17 18]);
model.component('comp1').physics('spf').create('inl1', 'InletBoundary', 1);
model.component('comp1').physics('spf').feature('inl1').selection.set([8]);
model.component('comp1').physics('spf').feature('inl1').set('BoundaryCondition', 'FullyDevelopedFlow');
model.component('comp1').physics('spf').feature('inl1').set('Uavfdf', 0.01);

model.component('comp1').multiphysics('tpf1').set('Fluid1', 'mat1');
model.component('comp1').multiphysics('tpf1').set('Fluid2', 'mat2');

model.component('comp1').physics('pf').create('iww1', 'InteriorWettedWall', 1);
model.component('comp1').physics('pf').feature('iww1').selection.set([4 5 9 10 14 16 17 18]);
model.component('comp1').physics('pf').feature('iww1').set('thetaw', '170[deg]');
model.component('comp1').physics('pf').create('inl1', 'InletBoundary', 1);
model.component('comp1').physics('pf').create('out1', 'Outlet', 1);
model.component('comp1').physics('pf').feature('out1').selection.set([1 3 20 21]);
model.component('comp1').physics('pf').feature('inl1').selection.set([8]);
model.component('comp1').physics('pf').feature('inl1').set('pfcond', 'Fluid2pf');

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 4);
model.component('comp1').mesh('mesh1').run('size');
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').selection.set([1 7]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 3);
model.component('comp1').mesh('mesh1').run('ftri1');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').run('ftri1');

model.study('std1').feature('time').set('tlist', 'range(0,0.005,2)');
model.study('std1').createAutoSequences('all');

model.result.dataset('dset1').set('geom', 'geom1');
model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').label([native2unicode(hex2dec({'90' '1f'}), 'unicode')  native2unicode(hex2dec({'5e' 'a6'}), 'unicode') ' (spf)']);
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').label([native2unicode(hex2dec({'88' '68'}), 'unicode')  native2unicode(hex2dec({'97' '62'}), 'unicode') ]);
model.result('pg1').feature('surf1').set('showsolutionparams', 'on');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('showsolutionparams', 'on');
model.result('pg1').feature('surf1').set('data', 'parent');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg2').label([native2unicode(hex2dec({'53' '8b'}), 'unicode')  native2unicode(hex2dec({'52' '9b'}), 'unicode') ' (spf)']);
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature.create('con1', 'Contour');
model.result('pg2').feature('con1').label([native2unicode(hex2dec({'7b' '49'}), 'unicode')  native2unicode(hex2dec({'50' '3c'}), 'unicode')  native2unicode(hex2dec({'7e' 'bf'}), 'unicode') ]);
model.result('pg2').feature('con1').set('showsolutionparams', 'on');
model.result('pg2').feature('con1').set('expr', 'p');
model.result('pg2').feature('con1').set('number', 40);
model.result('pg2').feature('con1').set('levelrounding', false);
model.result('pg2').feature('con1').set('smooth', 'internal');
model.result('pg2').feature('con1').set('showsolutionparams', 'on');
model.result('pg2').feature('con1').set('data', 'parent');
model.result.create('pg3', 'PlotGroup2D');
model.result('pg3').label([native2unicode(hex2dec({'6d' '41'}), 'unicode')  native2unicode(hex2dec({'4f' '53'}), 'unicode') ' 1 ' native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'4f' '53'}), 'unicode')  native2unicode(hex2dec({'79' 'ef'}), 'unicode')  native2unicode(hex2dec({'52' '06'}), 'unicode')  native2unicode(hex2dec({'65' '70'}), 'unicode') ' (pf)']);
model.result('pg3').set('frametype', 'spatial');
model.result('pg3').feature.create('surf1', 'Surface');
model.result('pg3').feature('surf1').set('expr', 'pf.Vf1');
model.result('pg3').feature('surf1').set('smooth', 'internal');
model.result('pg3').feature('surf1').set('data', 'parent');
model.result('pg3').feature.create('con1', 'Contour');
model.result('pg3').feature('con1').set('expr', 'pf.Vf1');
model.result('pg3').feature('con1').set('levelmethod', 'levels');
model.result('pg3').feature('con1').set('levels', '0.5');
model.result('pg3').feature('con1').set('coloring', 'uniform');
model.result('pg3').feature('con1').set('colorlegend', false);
model.result('pg3').feature('con1').set('color', 'gray');
model.result('pg3').feature('con1').set('smooth', 'none');
model.result('pg3').feature('con1').set('data', 'parent');
model.result.create('pg4', 'PlotGroup2D');
model.result('pg4').set('data', 'dset1');
model.result('pg4').label([native2unicode(hex2dec({'52' 'a8'}), 'unicode')  native2unicode(hex2dec({'7f' '51'}), 'unicode')  native2unicode(hex2dec({'68' '3c'}), 'unicode') ]);
model.result('pg4').create('mesh1', 'Mesh');
model.result('pg4').feature('mesh1').set('meshdomain', 'surface');
model.result('pg4').feature('mesh1').set('colortable', 'TrafficFlow');
model.result('pg4').feature('mesh1').set('colortabletrans', 'nonlinear');
model.result('pg4').feature('mesh1').set('nonlinearcolortablerev', true);
model.result('pg4').feature('mesh1').create('sel1', 'MeshSelection');
model.result('pg4').feature('mesh1').feature('sel1').selection.set([1 2 3 4]);
model.result('pg4').feature('mesh1').set('qualmeasure', 'custom');
model.result('pg4').feature('mesh1').set('qualexpr', 'comp1.spatial.relVol');
model.result('pg4').feature('mesh1').set('colorrangeunitinterval', false);
model.result.remove('pg2');
model.result.remove('pg1');
model.result.remove('pg4');
model.result.remove('pg3');

model.study('std1').createAutoSequences('sol');
model.study('std1').createAutoSequences('jobs');

model.result.dataset('dset1').set('geom', 'geom1');
model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').label([native2unicode(hex2dec({'90' '1f'}), 'unicode')  native2unicode(hex2dec({'5e' 'a6'}), 'unicode') ' (spf)']);
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').label([native2unicode(hex2dec({'88' '68'}), 'unicode')  native2unicode(hex2dec({'97' '62'}), 'unicode') ]);
model.result('pg1').feature('surf1').set('showsolutionparams', 'on');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('showsolutionparams', 'on');
model.result('pg1').feature('surf1').set('data', 'parent');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg2').label([native2unicode(hex2dec({'53' '8b'}), 'unicode')  native2unicode(hex2dec({'52' '9b'}), 'unicode') ' (spf)']);
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature.create('con1', 'Contour');
model.result('pg2').feature('con1').label([native2unicode(hex2dec({'7b' '49'}), 'unicode')  native2unicode(hex2dec({'50' '3c'}), 'unicode')  native2unicode(hex2dec({'7e' 'bf'}), 'unicode') ]);
model.result('pg2').feature('con1').set('showsolutionparams', 'on');
model.result('pg2').feature('con1').set('expr', 'p');
model.result('pg2').feature('con1').set('number', 40);
model.result('pg2').feature('con1').set('levelrounding', false);
model.result('pg2').feature('con1').set('smooth', 'internal');
model.result('pg2').feature('con1').set('showsolutionparams', 'on');
model.result('pg2').feature('con1').set('data', 'parent');
model.result.create('pg3', 'PlotGroup2D');
model.result('pg3').label([native2unicode(hex2dec({'6d' '41'}), 'unicode')  native2unicode(hex2dec({'4f' '53'}), 'unicode') ' 1 ' native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'4f' '53'}), 'unicode')  native2unicode(hex2dec({'79' 'ef'}), 'unicode')  native2unicode(hex2dec({'52' '06'}), 'unicode')  native2unicode(hex2dec({'65' '70'}), 'unicode') ' (pf)']);
model.result('pg3').set('frametype', 'spatial');
model.result('pg3').feature.create('surf1', 'Surface');
model.result('pg3').feature('surf1').set('expr', 'pf.Vf1');
model.result('pg3').feature('surf1').set('smooth', 'internal');
model.result('pg3').feature('surf1').set('data', 'parent');
model.result('pg3').feature.create('con1', 'Contour');
model.result('pg3').feature('con1').set('expr', 'pf.Vf1');
model.result('pg3').feature('con1').set('levelmethod', 'levels');
model.result('pg3').feature('con1').set('levels', '0.5');
model.result('pg3').feature('con1').set('coloring', 'uniform');
model.result('pg3').feature('con1').set('colorlegend', false);
model.result('pg3').feature('con1').set('color', 'gray');
model.result('pg3').feature('con1').set('smooth', 'none');
model.result('pg3').feature('con1').set('data', 'parent');
model.result.create('pg4', 'PlotGroup2D');
model.result('pg4').set('data', 'dset1');
model.result('pg4').label([native2unicode(hex2dec({'52' 'a8'}), 'unicode')  native2unicode(hex2dec({'7f' '51'}), 'unicode')  native2unicode(hex2dec({'68' '3c'}), 'unicode') ]);
model.result('pg4').create('mesh1', 'Mesh');
model.result('pg4').feature('mesh1').set('meshdomain', 'surface');
model.result('pg4').feature('mesh1').set('colortable', 'TrafficFlow');
model.result('pg4').feature('mesh1').set('colortabletrans', 'nonlinear');
model.result('pg4').feature('mesh1').set('nonlinearcolortablerev', true);
model.result('pg4').feature('mesh1').create('sel1', 'MeshSelection');
model.result('pg4').feature('mesh1').feature('sel1').selection.set([1 2 3 4]);
model.result('pg4').feature('mesh1').set('qualmeasure', 'custom');
model.result('pg4').feature('mesh1').set('qualexpr', 'comp1.spatial.relVol');
model.result('pg4').feature('mesh1').set('colorrangeunitinterval', false);
model.result.remove('pg2');
model.result.remove('pg1');
model.result.remove('pg4');
model.result.remove('pg3');

model.study('std1').feature('phasei').set('useadvanceddisable', true);
model.study('std1').feature('phasei').setSolveFor('/common/pres1', false);
model.study('std1').feature('phasei').set('disabledcommon', {'pres1'});
model.study('std1').createAutoSequences('sol');
model.study('std1').createAutoSequences('jobs');

model.sol('sol1').runFromTo('st1', 'v1');

model.result.dataset('dset1').set('geom', 'geom1');
model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').label([native2unicode(hex2dec({'90' '1f'}), 'unicode')  native2unicode(hex2dec({'5e' 'a6'}), 'unicode') ' (spf)']);
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').label([native2unicode(hex2dec({'88' '68'}), 'unicode')  native2unicode(hex2dec({'97' '62'}), 'unicode') ]);
model.result('pg1').feature('surf1').set('showsolutionparams', 'on');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('showsolutionparams', 'on');
model.result('pg1').feature('surf1').set('data', 'parent');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg2').label([native2unicode(hex2dec({'53' '8b'}), 'unicode')  native2unicode(hex2dec({'52' '9b'}), 'unicode') ' (spf)']);
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature.create('con1', 'Contour');
model.result('pg2').feature('con1').label([native2unicode(hex2dec({'7b' '49'}), 'unicode')  native2unicode(hex2dec({'50' '3c'}), 'unicode')  native2unicode(hex2dec({'7e' 'bf'}), 'unicode') ]);
model.result('pg2').feature('con1').set('showsolutionparams', 'on');
model.result('pg2').feature('con1').set('expr', 'p');
model.result('pg2').feature('con1').set('number', 40);
model.result('pg2').feature('con1').set('levelrounding', false);
model.result('pg2').feature('con1').set('smooth', 'internal');
model.result('pg2').feature('con1').set('showsolutionparams', 'on');
model.result('pg2').feature('con1').set('data', 'parent');
model.result.create('pg3', 'PlotGroup2D');
model.result('pg3').label([native2unicode(hex2dec({'6d' '41'}), 'unicode')  native2unicode(hex2dec({'4f' '53'}), 'unicode') ' 1 ' native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'4f' '53'}), 'unicode')  native2unicode(hex2dec({'79' 'ef'}), 'unicode')  native2unicode(hex2dec({'52' '06'}), 'unicode')  native2unicode(hex2dec({'65' '70'}), 'unicode') ' (pf)']);
model.result('pg3').set('frametype', 'spatial');
model.result('pg3').feature.create('surf1', 'Surface');
model.result('pg3').feature('surf1').set('expr', 'pf.Vf1');
model.result('pg3').feature('surf1').set('smooth', 'internal');
model.result('pg3').feature('surf1').set('data', 'parent');
model.result('pg3').feature.create('con1', 'Contour');
model.result('pg3').feature('con1').set('expr', 'pf.Vf1');
model.result('pg3').feature('con1').set('levelmethod', 'levels');
model.result('pg3').feature('con1').set('levels', '0.5');
model.result('pg3').feature('con1').set('coloring', 'uniform');
model.result('pg3').feature('con1').set('colorlegend', false);
model.result('pg3').feature('con1').set('color', 'gray');
model.result('pg3').feature('con1').set('smooth', 'none');
model.result('pg3').feature('con1').set('data', 'parent');
model.result.create('pg4', 'PlotGroup2D');
model.result('pg4').set('data', 'dset1');
model.result('pg4').setIndex('looplevel', 1, 0);
model.result('pg4').label([native2unicode(hex2dec({'52' 'a8'}), 'unicode')  native2unicode(hex2dec({'7f' '51'}), 'unicode')  native2unicode(hex2dec({'68' '3c'}), 'unicode') ]);
model.result('pg4').create('mesh1', 'Mesh');
model.result('pg4').feature('mesh1').set('meshdomain', 'surface');
model.result('pg4').feature('mesh1').set('colortable', 'TrafficFlow');
model.result('pg4').feature('mesh1').set('colortabletrans', 'nonlinear');
model.result('pg4').feature('mesh1').set('nonlinearcolortablerev', true);
model.result('pg4').feature('mesh1').create('sel1', 'MeshSelection');
model.result('pg4').feature('mesh1').feature('sel1').selection.set(+[1]);
model.result('pg4').feature('mesh1').set('qualmeasure', 'custom');
model.result('pg4').feature('mesh1').set('qualexpr', 'comp1.spatial.relVol');
model.result('pg4').feature('mesh1').set('colorrangeunitinterval', false);
model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').run;
model.result('pg3').feature('surf1').set('coloring', 'gradient');
model.result('pg3').run;
model.result('pg3').run;

model.sol('sol1').feature('t1').set('initialstepbdfactive', true);

model.study('std1').feature('time').set('plot', true);
model.study('std1').feature('time').set('plotfreq', 'tsteps');
model.study('std1').createAutoSequences('all');

model.component('comp1').physics('pf').feature('pfm1').set('epsilon_pf', '1.5*50[mm]/330');
model.component('comp1').physics('pf').feature('pfm1').set('chiOption', 'velocity');
model.component('comp1').physics('pf').feature('pfm1').set('U', '0.06[m/s]');

model.study('std1').setGenPlots(false);
model.study('std1').createAutoSequences('all');
model.study('std1').createAutoSequences('all');
model.study('std1').createAutoSequences('sol');
model.study('std1').createAutoSequences('jobs');

model.sol('sol1').runFromTo('st1', 'v1');

model.result('pg1').run;
model.result('pg3').run;
model.result('pg3').feature('surf1').set('coloring', 'gradient');
model.result('pg3').feature('surf1').set('colorlegend', false);

model.study('std1').feature('time').set('plotgroup', 'pg3');

model.component('comp1').physics('pf').feature('pfm1').set('U', 0.06);

model.study('std1').createAutoSequences('all');
model.study('std1').setGenPlots(true);
model.study('std1').createAutoSequences('all');

model.component('comp1').geom('geom1').feature.move('r5', 5);
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').run('mce1');
model.component('comp1').geom('geom1').run('uni2');

model.component('comp1').material('mat1').selection.set([1 2]);

model.label('Untitled12.mph');

model.component('comp1').common('pres1').set('analysisCaseEdited', true);
model.component('comp1').common('pres1').active(false);

model.component('comp1').physics('spf').active(false);
model.component('comp1').physics('spf').active(true);
model.component('comp1').physics('pf').feature('initfluid2').selection.set([2 3 4]);

model.component('comp1').multiphysics('tpf1').set('Mathematics_physics', 'none');
model.component('comp1').multiphysics('tpf1').set('IncludeSurfaceTension', false);

model.component('comp1').mesh('mesh1').run('ftri1');

model.study('std1').setGenPlots(false);
model.study('std1').createAutoSequences('all');

model.component('comp1').multiphysics('tpf1').set('Constitutiverelation1', 'Newtonian');

model.study('std1').createAutoSequences('all');

model.component('comp1').multiphysics('tpf1').set('IncludeSurfaceTension', false);
model.component('comp1').multiphysics('tpf1').set('Mathematics_physics', 'pf');

model.study('std1').createAutoSequences('all');

model.sol('sol1').runAll;

model.result('pg1').run;

model.component('comp1').geom('geom1').run('r5');
model.component('comp1').geom('geom1').run('r4');
model.component('comp1').geom('geom1').run('r4');
model.component('comp1').geom('geom1').run('uni2');

model.result('pg1').run;

model.label('Untitled12.mph');

model.component('comp1').common('pres1').set('analysisCaseEdited', true);

model.result('pg1').run;
model.result('pg2').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg3').run;
model.result('pg4').run;

model.component('comp1').geom('geom1').scaleUnitValue(true);
model.component('comp1').geom('geom1').scaleUnitValue(false);

model.component('comp1').physics('spf').feature('fp1').set('minput_temperature_src', 'fromCommonDef');
model.component('comp1').physics('spf').feature('open1').set('BoundaryCondition', 'NormalStress');
model.component('comp1').physics('pf').feature('initfluid2').set('InitialValuesOption', 'SpecifyPhase');

model.component('comp1').multiphysics('tpf1').set('minput_temperature_src', 'userdef');
model.component('comp1').multiphysics('tpf1').set('Constitutiverelation1', 'Newtonian');
model.component('comp1').multiphysics('tpf1').set('Fluid1', 'mat1');
model.component('comp1').multiphysics('tpf1').set('Fluid2', 'mat2');
model.component('comp1').multiphysics('tpf1').set('mu2_mat', 'userdef');
model.component('comp1').multiphysics('tpf1').set('mu1_mat', 'userdef');
model.component('comp1').multiphysics('tpf1').set('SurfaceTensionCoefficient', 'LibraryCoefficientLiquidGas');
model.component('comp1').multiphysics('tpf1').set('LibraryCoefficientLiquidGas', 'WaterAir');

model.study('std1').createAutoSequences('all');

app.form('form4').formObject('collection1').setIndex('vertscroll', false, 0);
app.form('form4').formObject('collection1').setIndex('vertscroll', true, 0);

model.label([native2unicode(hex2dec({'6c' '5f'}), 'unicode')  native2unicode(hex2dec({'97' '52'}), 'unicode')  native2unicode(hex2dec({'68' 'a6'}), 'unicode')  native2unicode(hex2dec({'67' '08'}), 'unicode') '.mph']);

model.component('comp1').common('pres1').set('analysisCaseEdited', true);

model.label([native2unicode(hex2dec({'6c' '5f'}), 'unicode')  native2unicode(hex2dec({'97' '52'}), 'unicode')  native2unicode(hex2dec({'68' 'a6'}), 'unicode')  native2unicode(hex2dec({'67' '08'}), 'unicode') '.mph']);

model.component('comp1').common('pres1').set('analysisCaseEdited', true);

out = model;

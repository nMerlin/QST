function tests = coherentStateTest
% COHERENTSTATETEST Verifies the correct behavior of COHERENTSTATE
tests = functiontests(localfunctions);
end

function testNormalization(testCase)
% Trace of density matrix should be 1
expSolution = 1;

% Fetch actual trace
rho = coherentState(30,5);
actSolution = trace(rho);

verifyEqual(testCase,actSolution,expSolution,'AbsTol',100*eps);
end
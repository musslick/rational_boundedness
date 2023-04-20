function [dualTaskCombs, similarity] = recomputeTaskDependencies(A_tasksIdx, dualTaskCombs, R_hidden, tasksToPerform)

            % test MSE of individual tasks depending on dependency type 3 (indirect, asummetric)
            
            % dependency types
            %
            %  O    O
            %  |      |
            %  O    O   type 0 (independent)
            %
            %
            %   O       type 1 (fan-in)
            %  / \   
            %
            %  \ /
            %   O       type 2 (fan-out)
            %
            %  O    O
            %  |  \  |
            %  O    O   type 3 (indirect, asymmetric)
            %
            %  O    O
            %  |  /\ |
            %  O    O   type 4 (indirect, symmetric)
            %
            
            similarity = nan(2, size(dualTaskCombs,2));
            
            % for each condition, check type of dependency
            
            for dualCondIdx = 1:size(dualTaskCombs,2)
                
                task1 = dualTaskCombs(1, dualCondIdx);
                task2 = dualTaskCombs(2, dualCondIdx);
                
                % determine dependency type
                [rowTask1, colTask1] = find(A_tasksIdx == task1);
                [rowTask2, colTask2] = find(A_tasksIdx == task2);
                
                if(colTask1 == colTask2 && rowTask1 ~= rowTask2) % fan-in dependency
                    dualTaskCombs(3, dualCondIdx) = 1;
                end
                
                if(rowTask1 == rowTask2 && colTask1 ~= colTask2) % fan-out dependency
                    dualTaskCombs(3, dualCondIdx) = 2;
                end
                
                if(rowTask1 ~= rowTask2 && colTask1 ~= colTask2 ...
                        && ((A_tasksIdx(rowTask2, colTask1) ~= 0 ...
                              && A_tasksIdx(rowTask1, colTask2) == 0) )) % asymmetric indirect dependency
                    dualTaskCombs(3, dualCondIdx) = 3;
                end
                
                if(rowTask1 ~= rowTask2 && colTask1 ~= colTask2 ...
                        && A_tasksIdx(rowTask2, colTask1) == 0 ...
                        && A_tasksIdx(rowTask1, colTask2) ~= 0) % asymmetric indirect dependency (mirrored)
                    dualTaskCombs(3, dualCondIdx) = 3;
                    dualTaskCombs(1:2, dualCondIdx) = flipud(dualTaskCombs(1:2, dualCondIdx));
                end
                
                if(rowTask1 ~= rowTask2 && colTask1 ~= colTask2 ...
                        && A_tasksIdx(rowTask2, colTask1) ~= 0 ...
                        && A_tasksIdx(rowTask1, colTask2) ~= 0) % symmetric indirect dependency
                    dualTaskCombs(3, dualCondIdx) = 4;
                end
                
                if(rowTask1 ~= rowTask2 && colTask1 ~= colTask2 ...
                        && A_tasksIdx(rowTask2, colTask1) == 0 ...
                        && A_tasksIdx(rowTask1, colTask2) == 0) % independent tasks
                    dualTaskCombs(3, dualCondIdx) = 0;
                    dualTaskCombs(1:2, dualCondIdx) = dualTaskCombs(randperm(2), dualCondIdx);  % randomly assign task 1 & task 2
                end
                
                
                
                % compute similairty
                task1 = dualTaskCombs(1, dualCondIdx);
                task2 = dualTaskCombs(2, dualCondIdx);
                
                [rowTask1, colTask1] = find(A_tasksIdx == task1);
                [rowTask2, colTask2] = find(A_tasksIdx == task2);
                
                sharedWithTask1 = A_tasksIdx(rowTask1, colTask2);
                sharedWithTask2 = A_tasksIdx(rowTask2, colTask1);
                
                if(sharedWithTask1 ~= 0)
                    indexTask1 = find(tasksToPerform == task1);
                    indexSharedWithTask1 = find(tasksToPerform == sharedWithTask1);
                    similarity(1, dualCondIdx) = R_hidden(indexTask1, indexSharedWithTask1);
                end
                
                if(sharedWithTask2 ~= 0)
                    indexTask2 = find(tasksToPerform == task2);
                    indexSharedWithTask2 = find(tasksToPerform == sharedWithTask2);
                    similarity(2, dualCondIdx) = R_hidden(indexTask2, indexSharedWithTask2);
                end
                
            end

            similarity = transpose(similarity);
end
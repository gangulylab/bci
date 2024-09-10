function key = constrain1D(targetID, decoder_out)
%function to constrain the movement of the cursor to one direction - towards the current target. 

% switch targetID
%     case 1
%         if decoder_out == 2
%             key = 2;
%         else key = 0;
%         end
%     case 3
%         if decoder_out == 3
%             key = 3;
%         else key = 0;
%         end
%     case 5
%         if decoder_out == 4
%             key = 4;
%         else key = 0;
%         end
%     case 7
%         if decoder_out == 1
%             key = 1;
%         else key = 0;
%         end
% end
if decoder_out ~= 7

switch targetID
    case 1
        if decoder_out ~= 0
            key = 5;
        else; key = 0;
        end
    case 3
        if decoder_out ~= 0
            key = 3;
        else; key = 0;
        end
    case 5
        if decoder_out ~= 0
            key = 6;
        else; key = 0;
        end
    case 7
        if decoder_out ~= 0
            key = 1;
        else; key = 0;
        end
end
else; key = 7;
end
end

